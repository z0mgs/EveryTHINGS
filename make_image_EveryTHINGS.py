# This script concat the ms files.
# Perform ms transform, and make the final image and moment maps
import sys
import imp
# from shutil import rmtree
from casac import casac
pipepath = '/data/scratch/idchiang/EveryTHINGS/'
# Import 3rd-party scripts
eiu = imp.load_source('everythings_image_utils',
                      pipepath + 'everythings_image_utils.py')
# Avoid sys package syntax error
print('#\nExecuted with Python ' + sys.version + '\n#\n')

msmd = casac.msmetadata()

# pbmask / pblimit
pbs = {'IC342': [0.15, 0.15],
       'NGC2787': [0.6, 0.2],
       'NGC3227': [0.8, 0.2],
       'NGC3608': [0.5, 0.2],
       'NGC3898': [0.8, 0.2],
       'NGC3941': [0.8, 0.7],
       'NGC4365': [0.2, 0.2],
       'NGC4374': [0.9, 0.2],
       'NGC4450': [0.8, 0.7],
       'NGC4477': [0.5, 0.2],
       'NGC4494': [0.5, 0.2],
       'NGC4496A': [0.9, 0.2],
       'NGC4501': [0.85, 0.2],
       'NGC4596': [0.5, 0.2],
       'NGC4636': [0.5, 0.2],
       'NGC4649': [0.5, 0.2],
       'NGC5728': [0.5, 0.2]}

###############################################################################
#
# The pipeline starts here
#
###############################################################################

trial_name_end = ''

#
# Read the info of the galaxy and list of ms files
#
target, mss, final_vis, glx_ctr = eiu.target_info(trial_name_end)
trial_name = target + trial_name_end
# Grab restfreq
msmd.open(mss[0])
freqs = msmd.chanfreqs(spw=0)
restfreq_Hz_num = (freqs[0] + freqs[-1]) / 2
msmd.close()
restfreq = str(round(restfreq_Hz_num / 1000000, 4)) + 'MHz'
del freqs
print('#\nTarget galaxy: ' + target + '\n' +
      'restfreq: ' + restfreq + '\n' +
      'Phase center: ' + glx_ctr + '\n#\n')
print('#\nFinal vis name: ' + final_vis + '\n#\n')

#
# w and nchan
#
if target in ['IC342']:
    w = 2.5
    nchan = 128
elif target in ['NGC3898', 'NGC4038', 'NGC4501']:
    w = 2.5
    nchan = 256
elif target in ['NGC4496A', 'NGC3631', 'NGC4380']:
    w = 5.0
    nchan = 64
elif target in ['NGC1961', 'NGC3227', 'NGC7479']:
    w = 5.0
    nchan = 128
elif target in ['NGC2787']:
    w = 10.0
    nchan = 64
elif target in ['NGC4941']:
    w = 4.0
    nchan = 128
else:
    w = 5.0
    nchan = 128

#
# Calculate image dimensions
#
try:
    pbmask, pblimit = pbs[target]
except KeyError:
    pbmask, pblimit = 0.4, 0.2

magnify_ratio = min((1 - pblimit) * 1.2, 1)

cellsize_str, cellsize_arcsec, imsize = \
    eiu.image_dimensions(vis=final_vis,
                         oversamplerate=5,
                         magnify_ratio=magnify_ratio)
print('#\nCalculated cell size: ' + cellsize_str + '\n' +
      'Image size: [' + str(imsize[0]) + ', ' + str(imsize[1]) + ']\n')

#
# Begin imageing
#
scales = [0, 5, 20, 70]
if target != 'IC342':
    glx_ctr = ''

#
eiu.imaging(vis=final_vis,
            trial_name=trial_name,
            interactive=False,
            # Image dimension
            imsize=imsize,
            cellsize=cellsize_str,
            glx_ctr=glx_ctr,
            restfreq=restfreq,
            specmode='cube',
            outframe='LSRK',
            veltype='radio',
            # Restore to common beam?
            restoringbeam='common',
            # Weighting
            weighting='briggs',
            robust=0.5,
            # Methods
            scales=scales,
            gain=0.1,
            smallscalebias=0.6,
            dogrowprune=False,
            growiterations=0,
            noisethreshold=2.0,  # Same as Dyas
            minbeamfrac=3.0,  # Same as Dyas
            sidelobethreshold=1.0,
            gridder='mosaic',
            pbmask=pbmask,
            pblimit=pblimit,
            # Stopping criteria
            threshold='1.0mJy',
            niter_in=5000000,  # don't want to hit
            nsigma=2.5,  # 1.0 too deep, 3.0 with cyclefactor 1.0 too shallow
            cyclefactor=1.0,
            minpsffraction=0.2,
            w=w,
            nchan=nchan)
