"""
Script name: casa_pipeline_lines_EveryTHINGS.py
Description: Script designed for JVLA HI 21cm emission line data calibration
             with CASA 5.4.1
Authors: I-Da Chiang, Dyas Utomo and the z0mgs team
E-mail: idchiang@ucsd.edu
Date: Apr. 22, 2019.

Reference 1: VLA CASA Calibration Pipeline 5.4.1
The general purpose JVLA reduction pipeline.
https://science.nrao.edu/facilities/vla/data-processing/pipeline

Reference 2: Astroua's LocalGroup-VLA pipeline
We refer heavily to their method of customized bandpass calibration
https://github.com/Astroua/LocalGroup-VLA/blob/master/14A-235/
pipeline_scripts/casa_pipeline_lines.py

Note: Put this script in the directory with sdm file
Note: Will need to change "pipepath" variable here and in
      everythings_manual_flag.py
Note: To run the pipeline, start casa environment with $casa --pipeline
"""
# General imports
import os
import sys
import tarfile
import traceback
from shutil import rmtree
# CASA imports
from h_init_cli import h_init_cli as h_init
from hifv_importdata_cli import hifv_importdata_cli as hifv_importdata
from hifv_hanning_cli import hifv_hanning_cli as hifv_hanning
from hifv_flagdata_cli import hifv_flagdata_cli as hifv_flagdata
from hifv_vlasetjy_cli import hifv_vlasetjy_cli as hifv_vlasetjy
from hifv_priorcals_cli import hifv_priorcals_cli as hifv_priorcals
# from hif_refant_cli import hif_refant_cli as hif_refant
from hifv_testBPdcals_cli import hifv_testBPdcals_cli as hifv_testBPdcals
from hifv_flagbaddef_cli import hifv_flagbaddef_cli as hifv_flagbaddef
# from hifv_uncalspw_cli import hifv_uncalspw_cli as hifv_uncalspw
from hifv_checkflag_cli import hifv_checkflag_cli as hifv_checkflag
from hifv_semiFinalBPdcals_cli import hifv_semiFinalBPdcals_cli as \
    hifv_semiFinalBPdcals
from hifv_solint_cli import hifv_solint_cli as hifv_solint
from hifv_fluxboot_cli import hifv_fluxboot_cli as hifv_fluxboot
from hifv_finalcals_cli import hifv_finalcals_cli as hifv_finalcals
from hifv_circfeedpolcal_cli import hifv_circfeedpolcal_cli as \
    hifv_circfeedpolcal
from hifv_applycals_cli import hifv_applycals_cli as hifv_applycals
from hifv_targetflag_cli import hifv_targetflag_cli as hifv_targetflag
# from hifv_statwt_cli import hifv_statwt_cli as hifv_statwt
from hifv_plotsummary_cli import hifv_plotsummary_cli as hifv_plotsummary
from hif_makeimlist_cli import hif_makeimlist_cli as hif_makeimlist
from hif_makeimages_cli import hif_makeimages_cli as hif_makeimages
from hifv_exportdata_cli import hifv_exportdata_cli as hifv_exportdata
from h_save_cli import h_save_cli as h_save
import pipeline.infrastructure.casatools as casatools
from tasks import split, uvcontsub, vishead
# EveryTHINGS imports
pipepath = '/data/scratch/idchiang/EveryTHINGS'  # EveryTHINGS pipeline path
sys.path.append(pipepath)
from everythings_utils import everythings_log, everythings_mwhi, \
    everythings_flagcali, everythings_flagtarget, everythings_statwt
from everythings_manual_flag import everythings_manual_flag
from everythings_diagnostic_plots import everythings_diagnostic_plots
from everythings_bandpass import astroua_bandpass

"""
# Setup paths; Should no longer be needed
sys.path.insert(0, os.path.expandvars("$SCIPIPE_HEURISTICS"))
execfile(os.path.join(os.path.expandvars("$SCIPIPE_HEURISTICS"),
                      "pipeline/h/cli/h.py"))
execfile(os.path.join(os.path.expandvars("$SCIPIPE_HEURISTICS"),
                      "pipeline/hif/cli/hif.py"))
execfile(os.path.join(os.path.expandvars("$SCIPIPE_HEURISTICS"),
                      "pipeline/hifa/cli/hifa.py"))
execfile(os.path.join(os.path.expandvars("$SCIPIPE_HEURISTICS"),
                      "pipeline/hifv/cli/hifv.py"))
"""

# Ignore any as refants??
refantignore = ""

# hifv general variables
importonly = False
hanning = False
polarization_cali = False
contsub = True
pipelinemode = 'automatic'
interactive = True
echo_to_screen = interactive
IMPORT_ONLY = ''
__rethrow_casa_exceptions = True

# Remove files from previous run
os.system('rm -r *.last *.restore *.txt *.fluxdensities *.ms.tar.gz')
os.system('rm -r weblog.* *.g *.k *.ms *.split *.png pipeline*')
os.system('rm -r *.ms.flagversions *.csv *.b final_caltables')

# EveryTHINGS data management parameters
remove_sdm = True
tar_ms = False
bandpass_fillgaps_n = 202

# Run the procedure
casatools.post_to_log("Beginning VLA pipeline run ...")

# Read SDM name
tgzname = ''
SDM_name = ''

filelist = os.listdir(os.getcwd())
for fn in filelist:
    temp = fn.split('.')
    if len(temp) == 5:
        SDM_name = fn
    elif (len(temp) == 8) and (temp[-2:] == ['tar', 'gz']):
        tgzname = fn
        if SDM_name == '':
            SDM_name = '.'.join(temp[:5])

assert SDM_name != '', "::EveryTHINGS:: Can't define valid SDM name."

if not os.path.isdir(SDM_name):
    assert tgzname != '', "::EveryTHINGS:: Can't find valid SDM.tar.gz file!!"
    with tarfile.open(tgzname) as tar:
        everythings_log("Start untar " + tgzname)
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tar)
        everythings_log("Finished untar " + tgzname)

vislist = [SDM_name]
msname = SDM_name + '.ms'

try:
    #
    # Initialize the pipeline
    #
    context = h_init(plotlevel='summary')
    #
    # Load the data
    #
    hifv_importdata(vis=vislist, pipelinemode=pipelinemode)
    if importonly:
        raise Exception(IMPORT_ONLY)
    if remove_sdm:
        everythings_log("Removing the SDM file")
        rmtree(SDM_name)
        everythings_log("The SDM file removed")
    #
    # Hanning smooth the data
    #
    if hanning:
        hifv_hanning(pipelinemode=pipelinemode)
    #
    # Quack
    #
    # Might need to enhance the quack here, depending on the first run.
    # Check out idchiang_HI_flagquack.py when needed
    #
    # Flag known bad data
    #
    hifv_flagdata(pipelinemode=pipelinemode, scan=True, hm_tbuff='1.5int',
                  intents='*POINTING*,*FOCUS*,*ATMOSPHERE*,' +
                  '*SIDEBAND_RATIO*,*UNKNOWN*,*SYSTEM_CONFIGURATION*,' +
                  '*UNSPECIFIED#UNSPECIFIED*')
    #
    # (EveryTHINGS) Flag channels with MW HI
    #
    everythings_mwhi(msname=msname, spw_HI=8, restfreq=1420405752,
                     numChanPhase=200, numChanBandpass=100)
    #
    # (EveryTHINGS) Flag manually identified defects
    #
    everythings_manual_flag(SDM_name=SDM_name, msname=msname)
    #
    # Fill model columns for primary calibrators
    #
    hifv_vlasetjy(pipelinemode=pipelinemode)
    #
    # Gain curves, opacities, antenna position corrections,
    # requantizer gains (NB: requires CASA 4.1!)
    #
    hifv_priorcals(pipelinemode=pipelinemode)
    #
    # Initial test calibrations using bandpass and delay calibrators
    #
    hifv_testBPdcals(pipelinemode=pipelinemode, refantignore=refantignore)
    astroua_bandpass(myvis=msname, stage='test', context=context,
                     refantignore=refantignore)
    #
    # Identify and flag basebands with bad deformatters or rfi based on
    # bp table amps and phases
    #
    hifv_flagbaddef(pipelinemode=pipelinemode)
    #
    # Flag spws that have no calibration at this point
    #
    # hifv_uncalspw(pipelinemode=pipelinemode, delaycaltable='testdelay.k',
    #               bpcaltable='testBPcal.b')
    #
    # Flag possible RFI on BP calibrator using rflag
    #
    hifv_checkflag(pipelinemode=pipelinemode)
    #
    # DO SEMI-FINAL DELAY AND BANDPASS CALIBRATIONS
    # (semi-final because we have not yet determined the spectral index of
    # the bandpass calibrator)
    #
    hifv_semiFinalBPdcals(pipelinemode=pipelinemode, refantignore=refantignore)
    #
    # Use flagdata rflag mode again on calibrators
    #
    hifv_checkflag(pipelinemode=pipelinemode, checkflagmode='semi')
    #
    # Re-run semi-final delay and bandpass calibrations
    #
    hifv_semiFinalBPdcals(pipelinemode=pipelinemode, refantignore=refantignore)
    astroua_bandpass(myvis=msname, stage='semifinal', context=context,
                     refantignore=refantignore)
    everythings_flagcali(msname=msname, spw='8~12')
    #
    # Flag spws that have no calibration at this point
    # hifv_uncalspw(pipelinemode=pipelinemode, delaycaltable='delay.k',
    #               bpcaltable='BPcal.b')
    #
    # Determine solint for scan-average equivalent
    #
    hifv_solint(pipelinemode=pipelinemode, refantignore=refantignore)
    #
    # Do the flux density boostrapping -- fits spectral index of
    # calibrators with a power-law and puts fit in model column
    #
    hifv_fluxboot(pipelinemode=pipelinemode, refantignore=refantignore)
    #
    # Make the final calibration tables
    #
    hifv_finalcals(pipelinemode=pipelinemode, refantignore=refantignore)
    astroua_bandpass(myvis=msname, stage='final', context=context,
                     refantignore=refantignore)
    #
    # Polarization calibration
    #
    if polarization_cali:
        hifv_circfeedpolcal(pipelinemode=pipelinemode)
    #
    # Apply all the calibrations and check the calibrated data
    #
    hifv_applycals(pipelinemode=pipelinemode)
    #
    # Flag spws that have no calibration at this point
    #
    # hifv_uncalspw(pipelinemode=pipelinemode, delaycaltable='finaldelay.k',
    #               bpcaltable='finalBPcal.b')
    #
    # Now run all calibrated data, including the target, through rflag
    #
    hifv_targetflag(pipelinemode=pipelinemode, intents='*CALIBRATE*')
    everythings_flagtarget(msname=msname, spw='8~12',
                           cont_spw='8:0~1024;3072~4096')
    #
    # Calculate data weights based on standard deviation within each spw
    #
    # hifv_statwt(pipelinemode=pipelinemode)
    everythings_statwt(msname=msname, cont_spw='8:0~1024;3072~4096')
    #
    # Plotting Summary
    #
    hifv_plotsummary(pipelinemode=pipelinemode)
    #
    # Make a list of expected point source calibrators to be cleaned
    #
    hif_makeimlist(intent='PHASE,BANDPASS', specmode='cont', spw='8',
                   pipelinemode=pipelinemode)
    #
    # Make clean images for the selected calibrators
    #
    hif_makeimages(hm_masking='none')
    #
    # Export the data
    #
    os.mkdir('products/')
    hifv_exportdata(products_dir='products/')  # I don't know what will happen

except Exception as e:
    if str(e) == IMPORT_ONLY:
        casatools.post_to_log("Exiting after import step ...",
                              echo_to_screen=echo_to_screen)
    else:
        casatools.post_to_log("Error in procedure execution ...",
                              echo_to_screen=echo_to_screen)
        errstr = traceback.format_exc()
        casatools.post_to_log(errstr, echo_to_screen=echo_to_screen)

finally:
    #
    # Save the results to the context
    #
    h_save()
    #
    # My own diagnostic plots
    #
    # Copyfile not finished
    everythings_diagnostic_plots(msname=msname, plot_spw=8)
    #
    # uvcontsub
    #
    msforsplit = msname
    dataforsplit = 'corrected'
    spwsplit = '8'
    if contsub:
        everythings_log("Start uvcontsub")
        uvcontsub(vis=msname, fitspw='8:1024~3072', excludechans=True,
                  want_cont=False, fitorder=1, spw='8')
        msforsplit = msname + '.contsub'
        dataforsplit = 'data'
        spwsplit = ''
    #
    # Split
    #
    everythings_log("Start final split")
    mssplit = msname + '.split'
    split(vis=msforsplit, outputvis=mssplit, spw=spwsplit,
          intent='*TARGET*', width=4, datacolumn=dataforsplit)
    # Set spw name to avoid further warning
    temp = vishead(mssplit, mode='get', hdkey='spw_name')
    temp = list(temp[0][0])
    temp[-1] = '0'
    temp = ''.join(temp)
    vishead(mssplit, mode='put', hdkey='spw_name', hdvalue=temp)
    # Remove contsub file
    if contsub:
        rmtree(msforsplit)
    #
    # Tar original ms
    #
    if tar_ms:
        source_dir = msname + ".flagversions"
        if os.path.isdir(source_dir):
            everythings_log("Remove .flagversions")
            rmtree(source_dir)
        everythings_log("Tar .ms")
        source_dir = msname
        if os.path.isdir(source_dir):
            with tarfile.open(source_dir + ".tar.gz", "w:gz") as tar:
                tar.add(source_dir, arcname=os.path.basename(source_dir))
            everythings_log("Remove .ms")
            rmtree(source_dir)
    #
    casatools.post_to_log("VLA CASA Pipeline finished." +
                          "  Terminating procedure execution ...",
                          echo_to_screen=echo_to_screen)
    #
    # Restore previous state
    #
    # __rethrow_casa_exceptions = def_rethrow
