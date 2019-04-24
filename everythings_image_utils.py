# General imports
import os
import sys
import gzip
import shutil
import imp
import numpy as np
# import matplotlib as mpl
# import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy import ndimage
from shutil import rmtree
from astropy.io import fits
import astropy.units as u
from astropy.convolution import convolve
from astropy.coordinates import Angle
# from shutil import copytree
# CASA imports
from tasks import imstat, immoments, imsmooth, imhead, immath, exportfits, \
    tclean, concat, imcontsub, vishead
from casac import casac
ia = casac.image()
# au imports
aupath = '/home/idchiang/bin/au/'
sys.path.append(aupath)
import analysisUtils as au
gal_data = \
    imp.load_source('gal_data', '/home/idchiang/idc_lib/idc_lib/gal_data.py')


def target_info(trial_name_end=''):
    # Read and return target info
    target = os.getcwd().split('/')[-1]
    trial_name = target + trial_name_end
    os.system('rm -r ' + trial_name + '.*')
    subdirs = ['01', '02', '03']
    mss = []
    for subdir in subdirs:
        if os.path.isdir(subdir):
            fns = os.listdir(subdir)
            for fn in fns:
                temp = fn.split('.')
                if len(temp) > 4:
                    if temp[-1] == 'split':
                        mss.append(subdir + '/' + fn)
    # Calculate restfreq and phasecenter
    gdata = gal_data.gal_data(target)
    ra = str(round(gdata.field('RA_DEG')[0], 1)) + 'deg '
    dec = str(round(gdata.field('DEC_DEG')[0], 1)) + 'deg'
    glx_ctr = 'J2000 ' + ra + dec
    # target: name of target galaxy
    # mss: the calibrated ms files
    final_vis = target + '.ms'
    if os.path.isdir(final_vis):
        print("...Removing previous final vis file")
        rmtree(final_vis)
    print("...Concating " + str(len(mss)) + " ms files for " +
          target)
    concat(vis=mss, concatvis=final_vis)
    return target, mss, final_vis, glx_ctr


def image_dimensions(vis, oversamplerate=5, magnify_ratio=1.2):
    # Dimensions from analysisUtils
    au_cellsize, au_imsize, au_oversample = \
        au.pickCellSize(vis=vis, imsize=True, npix=oversamplerate)
    # Determine cellsize
    if au_cellsize < 1.0:
        cellsize_arcsec = 1.0
    else:
        cellsize_arcsec = int(au_cellsize * 2) / 2.0
    cellsize_str = str(cellsize_arcsec) + 'arcsec'
    # Determin image size
    valid_sizes = []
    for ii in range(10):
        for kk in range(3):
            for jj in range(3):
                valid_sizes.append(2**(ii+1)*5**(jj)*3**(kk))
    valid_sizes.sort()
    valid_sizes = np.array(valid_sizes)
    need_cells_x = au_imsize[0] * magnify_ratio * au_cellsize / cellsize_arcsec
    need_cells_y = au_imsize[1] * magnify_ratio * au_cellsize / cellsize_arcsec
    cells_x = np.min(valid_sizes[valid_sizes > need_cells_x])
    cells_y = np.min(valid_sizes[valid_sizes > need_cells_y])
    imsize = [cells_x, cells_y]
    return cellsize_str, cellsize_arcsec, imsize


def make_mask_python(image_name, masked_beam=2.0, masked_spec_width=3,
                     masked_snr=2.0, reject_area_in_beams=2):
    """
    image_name: trial_name + '.image'
    masked_beam: ratio of convolved/origin bmaj/bmin
    masked_spec_width: num of channels to boxcar. Change in parent script
    masked_snr: sigma level to keep/remove
    reject_area_in_beams: Masked area rejected in unit of beam
    """
    image_spat_smo = image_name + '.spat_smo'
    image_spat_spec_smo = image_name + '.spat_spec_smo'
    image_snr_cut = image_name + '.spat_spec_smo.snr_cut'
    image_mask = image_name + '.pymask'
    #
    FWHM_TO_AREA = 2 * np.pi / (8 * np.log(2))
    # Adapted from Dyas and Adam's script
    print(" ::z0mgs:: Make mask for the image cube")
    #
    # 1) Spatially smooth the cube
    #
    # 1-1) Read bmaj, bmin, bpa
    temp_beam = imhead(image_name)['restoringbeam']
    bmaj = temp_beam['major']
    bmin = temp_beam['minor']
    bpa = temp_beam['positionangle']
    bmaj_smo = str(masked_beam*bmaj['value']) + bmaj['unit']
    bmin_smo = str(masked_beam*bmin['value']) + bmin['unit']
    bpa_smo = str(bpa['value']) + bpa['unit']
    # 1-2) Calculate beam area for later use
    assert imhead(image_name)['axisnames'][0] == 'Right Ascension'
    assert imhead(image_name)['axisnames'][1] == 'Declination'
    assert imhead(image_name)['axisunits'][0] == 'rad'
    assert imhead(image_name)['axisunits'][1] == 'rad'
    assert bmaj['unit'] == 'arcsec'
    assert bmin['unit'] == 'arcsec'
    incr_arcsec = Angle(np.mean(np.abs(imhead(image_name)['incr'][:2])) *
                        u.rad).arcsec
    bmaj_pix = bmaj['value'] / incr_arcsec
    bmin_pix = bmin['value'] / incr_arcsec
    beam_area_pix = bmaj_pix * bmin_pix * FWHM_TO_AREA
    # 1-3) Convolve the image
    print(' ::z0mgs:: Smoothed bmaj: ' + bmaj_smo)
    print(' ::z0mgs:: Smoothed bmin: ' + bmin_smo)
    print(' ::z0mgs:: Smoothed bpa: ' + bpa_smo)
    if os.path.isdir(image_spat_smo):
        rmtree(image_spat_smo)
    imsmooth(
        imagename=image_name,
        targetres=True,
        major=bmaj_smo,
        minor=bmin_smo,
        pa=bpa_smo,
        outfile=image_spat_smo,
        overwrite=True)
    #
    # 2) Spectrally smooth the cube
    #
    # 2-1) Locate the spectral axis
    im_head = imhead(image_spat_smo)
    try:
        freq_axis = np.argwhere(im_head['axisnames'] == 'Frequency')[0, 0]
    except IndexError:
        print(' ::z0mgs:: Cannot find the frequency axis!!')
        print(' ::z0mgs:: Here are the axisnames:')
        print(' ::z0mgs::', im_head['axisnames'])
        sys.exit()
    # 2-2) Read the casa image into numpy array
    ia.open(image_spat_smo)
    array_spat_smo = ia.getchunk()
    # 2-3) boxcar the spectral image according to spectral axis found
    array_spat_smo_swap = np.moveaxis(array_spat_smo, freq_axis, -1)
    array_spat_spec_smo_swap = np.empty_like(array_spat_smo_swap)
    s = array_spat_smo_swap.shape
    assert len(s) == 4  # Haven't design for other cases
    if masked_spec_width > 1:
        if masked_spec_width % 2:
            w = np.full(masked_spec_width, 1.0 / masked_spec_width)
        else:
            w = np.ones(masked_spec_width + 1, dtype=float)
            w[0] = 0.5
            w[-1] = 0.5
            w = w / w.sum()
        for i in range(s[0]):
            for j in range(s[1]):
                for k in range(s[2]):
                    array_spat_spec_smo_swap[i, j, k] = \
                        convolve(array=array_spat_smo_swap[i, j, k],
                                 kernel=w,
                                 boundary='extend')
    else:
        array_spat_spec_smo_swap = array_spat_smo_swap
    array_spat_spec_smo = np.moveaxis(array_spat_spec_smo_swap, -1, freq_axis)
    del array_spat_spec_smo_swap, array_spat_smo_swap
    if os.path.isdir(image_spat_spec_smo):
        rmtree(image_spat_spec_smo)
    ia2 = ia.subimage(outfile=image_spat_spec_smo)
    ia2.putchunk(array_spat_spec_smo)
    ia.close()
    ia2.close()
    #
    # 3) Estimate the noise level
    #
    # 3-1) For comparison, print the directly calculated rms values
    temp_stats = imstat(image_spat_smo)
    temp_rms = temp_stats['medabsdevmed'][0] / 0.6745
    print(' ::z0mgs:: RMS from spat_smo: ' + str(temp_rms))
    temp_stats = imstat(image_spat_spec_smo, algorithm='classic')
    temp_rms = temp_stats['medabsdevmed'][0] / 0.6745
    print(' ::z0mgs:: RMS from spat_spec_smo: ' + str(temp_rms))
    #
    # *) Remove these things first. Make sure I don't use them again
    #
    rmtree(image_spat_smo)
    del image_spat_smo, temp_stats, temp_rms
    #
    work_stats = imstat(image_spat_spec_smo, algorithm='chauvenet',
                        zscore=5)
    work_rms = work_stats['medabsdevmed'][0] / 0.6745
    print(' ::z0mgs:: RMS from spat_spec_smo: ' + str(work_rms))
    #
    # 3-2) Finding pixels above the noise level
    #
    if os.path.isdir(image_snr_cut):
        rmtree(image_snr_cut)
    immath(
        imagename=image_spat_spec_smo,
        outfile=image_snr_cut,
        expr='iif(IM0 >= ' + str(masked_snr*work_rms) + ',1.0,0.0)')
    #
    # 4) Reject small region
    #
    ia.open(image_snr_cut)
    mask = ia.getchunk()
    # ia.done()
    regions, n_regions = ndimage.label(mask)
    myhistogram = ndimage.measurements.histogram(regions, 0, n_regions + 1,
                                                 n_regions + 1)
    object_slices = ndimage.find_objects(regions)
    for i in range(n_regions):
        if myhistogram[i + 1] < reject_area_in_beams * beam_area_pix:
            mask[object_slices[i]] = 0
    #
    # 4-1) Convert the numpy array into casa image
    #
    if os.path.isdir(image_mask):
        rmtree(image_mask)
    ia2 = ia.subimage(outfile=image_mask)
    ia2.putchunk(mask)
    ia.close()
    ia2.close()
    #
    # Removing junks
    #
    rmtree(image_spat_spec_smo)
    rmtree(image_snr_cut)
    return image_mask, work_rms


def no_detection_check(image, pbimage, target, work_rms, chanwidth):
    # 1) Examine rms
    # 1-1) Convert the rms to corresponding noise level in Kelvin
    ia.open(image)
    data = ia.getchunk()
    ia.close()
    temp_beam = imhead(image)['restoringbeam']
    bmaj = temp_beam['major']
    bmin = temp_beam['minor']
    assert bmaj['unit'] == 'arcsec'
    assert bmin['unit'] == 'arcsec'
    print('::EveryTHINGS:: BMAJ: ' + str(bmaj['value']))
    print('::EveryTHINGS:: BMIN: ' + str(bmin['value']))
    noise_kelvin = 6.07 * 10**5 * work_rms / bmaj['value'] / bmin['value']
    noise_column = 1.823 * 10**18 * noise_kelvin * chanwidth
    #
    g = gal_data(target)
    try:
        morph = g['MORPH'][0]
    except KeyError:
        morph = 'morph not available'
    del g
    #
    work_rms = noise_column
    data *= 6.07 * 10**5 / bmaj['value'] / bmin['value']
    data *= 1.823 * 10**18 * chanwidth
    # work_rms *= 1000
    # data *= 1000
    #
    # Need: PB image
    # 2) Spectrum check
    # 2-0) Build PB mask
    ia.open(pbimage)
    pb = ia.getchunk()
    ia.close()
    mask99 = pb >= 0.99
    mask95 = pb >= 0.95
    mask90 = pb >= 0.90
    # 2-1) Plot the spectrum of the galaxy center
    spectrum99 = []
    spectrum95 = []
    spectrum90 = []
    numChan = data.shape[-1]
    for i in range(numChan):
        spectrum99.append(np.nanmean(data[..., i][mask99[..., i]]))
        spectrum95.append(np.nanmean(data[..., i][mask95[..., i]]))
        spectrum90.append(np.nanmean(data[..., i][mask90[..., i]]))
    spectrum99 = np.array(spectrum99)
    spectrum95 = np.array(spectrum95)
    spectrum90 = np.array(spectrum90)  # Not super meaningful
    fig, ax = plt.subplots()
    ax.plot(spectrum99, label='Spectrum-PB99')
    ax.plot(spectrum95, label='Spectrum-PB95')
    ax.plot(spectrum90, label='Spectrum-PB90')
    ax.plot([work_rms] * numChan, label=r'1$\sigma$')
    ax.plot([2 * work_rms] * numChan, label=r'2$\sigma$')
    ax.plot([3 * work_rms] * numChan, label=r'3$\sigma$')
    ax.set_ylabel(r"N$_{HI}$ [1/cm$^2$]")
    ax.set_xlabel("Channel")
    ax.set_title(target + ' (' + morph + ')')
    ax.legend()
    fig.savefig(image + 'noise_pb09spectrum.png')
    # 3) Output the results in a txt file
    """
    results = np.array([['Noise level:', str(noise_kelvin), 'Kelvin'],
                        ['Noise level:', '{:0.2e}'.format(noise_column),
                         '1/cm2']
                        ['SpectrumPB99 >1 sigma:',
                         str(np.sum(spectrum99 > work_rms)),
                         '/ ' + str(numChan)],
                        ['SpectrumPB99 >1 sigma:',
                         str(np.sum(spectrum99 > 2 * work_rms)),
                         '/ ' + str(numChan)],
                        ['SpectrumPB99 >1 sigma:',
                         str(np.sum(spectrum99 > 3 * work_rms)),
                         '/ ' + str(numChan)]])
    np.savetxt('noise_properties.txt', results, fmt='%s')
    """


def everythings_imcontsub(image='', fitorder=1, masked_spec_width=3,
                          pbimage='', w=0):
    assert False, 'Function under development!'
    # 0) Definitions
    imname_cont = image + '.imcont'
    imname_contsub = image + '.imcontsub'
    # 1) Load cube
    vh = vishead(image)
    numchan = vh  # Need to write 
    # 1-1) Select channels according to channel number
    chans = '0~10,' + str(numchan - 11) + '~' + str(numchan - 1)
    # 2) imcontsub. save both cont and contsub
    imcontsub(imagename=image,
              linefile=imname_contsub,
              contfile=imname_cont,
              fitorder=fitorder,
              chans=chans)
    # 3) run masking for contsub image with pbcor image itself
    mask_dir, work_rms = \
        make_mask_python(imname_contsub, masked_beam=2.0,
                         masked_spec_width=masked_spec_width,
                         masked_snr=2.0, reject_area_in_beams=2)
    # 4) run moment map building for contsub image
    cubename = image + '.imcs'
    tempcube = cubename + '.cube.fits'
    m0name = cubename + '.integrated'
    tempm0 = cubename + '.m0.fits'
    m1name = cubename + '.weighted_coord'
    tempm1 = cubename + '.m1.fits'
    m2name = cubename + '.weighted_dispersion_coord'
    tempm2 = cubename + '.m2.fits'
    for mxname in [m0name, m1name, m2name]:
        if os.path.isdir(mxname):
            rmtree(mxname)
    immoments(imagename=imname_contsub, moments=[0, 1, 2], axis='spectral',
              mask=mask_dir, excludepix=-1, outfile=cubename)
    for (mxname, fitsname) in zip([imname_contsub, m0name, m1name, m2name],
                                  [tempcube, tempm0, tempm1, tempm2]):
        exportfits(imagename=mxname, fitsimage=fitsname)
        # Blocking nan pixels
        if mxname == imname_contsub:
            data = fits.getdata(tempcube)
            nanmask = np.any(np.isnan(data), axis=(0, 1))
            del data
        if mxname in [m0name, m1name, m2name]:
            print("...Setting " + mxname + " NaNs to zeros")
            data, hdr = fits.getdata(fitsname, header=True)
            data[np.isnan(data)] = 0.0
            data[0, 0][nanmask] = np.nan
            fits.writeto(fitsname, data, hdr, overwrite=True)
        with open(fitsname, 'rb') as f_in:
            with gzip.open(fitsname + '.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(fitsname)
        if mxname in [m0name, m1name, m2name]:
            rmtree(mxname)
    # 5) run no detection check
    target = os.getcwd().split('/')[-1]
    no_detection_check(image=imname_contsub,
                       pbimage=pbimage,
                       target=target,
                       work_rms=work_rms,
                       chanwidth=w)


def imaging(vis, trial_name, interactive, imsize, cellsize, glx_ctr, restfreq,
            specmode, outframe, veltype, restoringbeam, weighting, robust,
            scales, smallscalebias, dogrowprune,
            growiterations, noisethreshold, minbeamfrac, sidelobethreshold,
            gridder, pbmask, pblimit, threshold, niter_in, nsigma, cyclefactor,
            minpsffraction, gain, w, nchan):
    #
    titles = ['.dirty', '', '.2.strong']
    #
    for i in range(2):
        # niter
        niter = 0 if i == 0 else niter_in
        # masking method
        usemask = 'auto-multithresh' if i == 2 else 'pb'
        # deconvolver
        deconvolver = 'multiscale'
        #
        tclean(vis=vis,
               imagename=trial_name,
               interactive=interactive,
               intent='*TARGET*',
               #
               datacolumn='data',
               nchan=nchan,
               start=str(w * nchan / (-2)) + 'km/s',
               width=str(w) + 'km/s',
               # Image dimension
               imsize=imsize,
               cell=cellsize,
               phasecenter=glx_ctr,
               restfreq=restfreq,
               specmode=specmode,
               outframe=outframe,
               veltype=veltype,
               # Restore to common beam?
               restoringbeam=restoringbeam,
               # Weighting
               weighting=weighting,
               robust=robust,
               # Methods
               deconvolver=deconvolver,
               scales=scales,
               gain=gain,
               smallscalebias=smallscalebias,
               usemask=usemask,
               dogrowprune=dogrowprune,
               growiterations=growiterations,
               noisethreshold=noisethreshold,
               minbeamfrac=minbeamfrac,
               sidelobethreshold=sidelobethreshold,
               gridder=gridder,
               pbmask=pbmask,
               pblimit=pblimit,
               pbcor=True,
               # Stopping criteria
               threshold=threshold,
               niter=niter,
               nsigma=nsigma,
               cyclefactor=cyclefactor,
               minpsffraction=minpsffraction)
        #
        nonpbcube = trial_name + '.image'
        pbcube = trial_name + '.image.pbcor'
        cubename = trial_name + '.image.pb'
        tempcube = trial_name + titles[i] + '.cube.fits'
        residualname = trial_name + '.residual'
        tempresidual = trial_name + titles[i] + '.residual.fits'
        m0name = cubename + '.integrated'
        tempm0 = trial_name + titles[i] + '.m0.fits'
        m1name = cubename + '.weighted_coord'
        tempm1 = trial_name + titles[i] + '.m1.fits'
        m2name = cubename + '.weighted_dispersion_coord'
        tempm2 = trial_name + titles[i] + '.m2.fits'
        for mxname in [m0name, m1name, m2name]:
            if os.path.isdir(mxname):
                rmtree(mxname)
        #
        print("...Calculating moment maps")
        if i == 0:
            rms = imstat(pbcube)['rms'][0]
            excludepix = [-1000, (2 * rms)]
            immoments(imagename=pbcube, moments=[0], axis='spectral',
                      excludepix=excludepix, outfile=m0name)
        else:
            # calculate mask here
            # change masked_spec_width according to w
            if w > 10:
                spec_width = 1
            else:
                spec_width = 3
            mask_dir, work_rms = \
                make_mask_python(nonpbcube, masked_beam=2.0,
                                 masked_spec_width=spec_width,
                                 masked_snr=2.0, reject_area_in_beams=2)
            # Use the PB corrected image!
            immoments(imagename=pbcube, moments=[0, 1, 2], axis='spectral',
                      mask=mask_dir, excludepix=-1, outfile=cubename)
        print("#\n")
        print("...Exporting fits files of " + titles[i])
        for (mxname, fitsname) in zip([pbcube, residualname, m0name, m1name,
                                       m2name],
                                      [tempcube, tempresidual, tempm0, tempm1,
                                       tempm2]):
            if (mxname in [residualname, m1name, m2name]) and i == 0:
                continue
            exportfits(imagename=mxname, fitsimage=fitsname)
            # Blocking nan pixels
            if mxname == pbcube:
                data = fits.getdata(tempcube)
                nanmask = np.any(np.isnan(data), axis=(0, 1))
                del data
            if mxname in [m0name, m1name, m2name]:
                print("...Setting " + mxname + " NaNs to zeros")
                data, hdr = fits.getdata(fitsname, header=True)
                data[np.isnan(data)] = 0.0
                data[0, 0][nanmask] = np.nan
                fits.writeto(fitsname, data, hdr, overwrite=True)
            with open(fitsname, 'rb') as f_in:
                with gzip.open(fitsname + '.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(fitsname)
            if mxname in [m0name, m1name, m2name]:
                rmtree(mxname)
        """
        print("...Plotting weighted m1 image")
        plt.ioff()
        m0 = fits.getdata(tempm0 + '.gz')[0, 0]
        m1 = fits.getdata(tempm1 + '.gz')[0, 0]
        #
        cmap = cm.bwr_r
        norm = mpl.colors.Normalize(vmin=np.nanmin(m1),
                                    vmax=np.nanmax(m1))
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        with np.errstate(invalid='ignore'):
            mm = m.to_rgba(m1)[:, :, :3]
            temp = np.sum(mm, axis=2)
            for i in range(3):
                mm[:, :, i] /= temp
            m0[m0 <= 0] = np.min(m0[m0 > 0])
        m0 = np.log10(m0)
        m0 -= m0.min()
        m0 /= m0.max()
        for i in range(3):
            mm[:, :, i] *= m0
        mm /= mm.max()
        #
        fig1, ax1 = plt.subplots(figsize=(12, 10))
        fig2, ax2 = plt.subplots()
        mpb = ax2.imshow(m1, cmap='bwr_r', origin='lower')
        ax1.imshow(mm, origin='lower')
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        plt.colorbar(mpb, ax=ax1)
        fig1.savefig(trial_name + titles[i] + '.m1.m0-weighted.png')
        plt.close('all')
        """
        if i == 1:
            target = os.getcwd().split('/')[-1]
            no_detection_check(image=pbcube,
                               pbimage=trial_name + '.pb',
                               target=target,
                               work_rms=work_rms,
                               chanwidth=w)
