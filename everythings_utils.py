"""
Authors: I-Da Chiang, Dyas Utomo and the z0mgs team
E-mail: idchiang@ucsd.edu
"""
# General imports
import numpy as np
# CASA imports
import pipeline.infrastructure.casatools as casatools
from casac import casac
from tasks import flagdata
from tasks import oldstatwt as statwt

# CASA tools
msmd = casac.msmetadata()


# Post message to casa log and terminal, with EveryTHINGS marker
def everythings_log(msg):
    casatools.post_to_log("::EveryTHINGS:: " + msg)


# Identify channels with MW HI
def everythings_mwhi(msname='', spw_HI=8, restfreq=1420405752,
                     numChanPhase=200, numChanBandpass=100):
    # spw_HI: int
    # restfreq: int. in Hz.
    # numChan: int
    msmd.open(msname)
    freqs = msmd.chanfreqs(spw_HI)
    freqs_delta = np.mean(freqs[1:] - freqs[:-1])
    # i=0: Phase calibrator
    # i=1: Bandpass calibrator:
    for i in range(2):
        numChan = numChanPhase if i == 0 else numChanBandpass
        min_freq = restfreq - freqs_delta * numChan
        max_freq = restfreq + freqs_delta * numChan
        mask = (freqs < max_freq) * (freqs > min_freq)
        if np.sum(mask) > 0:
            temp = np.arange(len(mask), dtype=int)
            chan_ids = temp[mask]
            min_chan, max_chan = \
                str(np.min(chan_ids)), str(np.max(chan_ids) + 1)
            rm_spw = str(spw_HI) + ':' + min_chan + '~' + max_chan
            intent = '*PHASE*' if i == 0 else '*BANDPASS*'
            everythings_log("Removing MW HI in " + intent + " calibrator.")
            flagdata(vis=msname, spw=rm_spw, intent=intent)
            del rm_spw, chan_ids, temp, min_chan, max_chan
        del mask, max_freq, min_freq
    msmd.close()


# Improved flagging for calibrators after semifinal BPcals
def everythings_flagcali(msname='', spw='8~12'):
    # Narrow-time, broad-band RFI flag at both calibrators at non-continuum
    # Extend flag is broken. So separate this to a two-stage flagging
    everythings_log("Extra flagging for calibrators")
    flagdata(vis=msname, intent='*CALI*', mode='tfcrop',
             datacolumn='corrected', spw=spw,
             extendflags=False, action='apply', display='',
             freqcutoff=3.0, timecutoff=2.0,
             combinescans=True,
             winsize=1, flagdimension='timefreq')
    flagdata(vis=msname, intent='*CALI*', mode='extend',
             datacolumn='corrected', spw=spw,
             extendflags=True, action='apply', display='',
             freqcutoff=100.0, timecutoff=100.0,
             growtime=70.0, growfreq=40.0,  # weak time due to quack
             extendpols=True, combinescans=False,
             winsize=1, flagdimension='timefreq')


# Target flag.
def everythings_flagtarget(msname='', spw='8~12',
                           cont_spw='8:0~1024;3072~4096'):
    everythings_log("Extra flagging for target")
    # time only
    msmd.open(msname)
    idc_target_scans = msmd.scansforintent('*TARGET*')
    idc_target_scans_str = [str(its) for its in idc_target_scans]
    msmd.close()
    for itss in idc_target_scans_str:
        flagdata(vis=msname, scan=itss, mode='tfcrop',
                 datacolumn='corrected', spw=spw,
                 extendflags=False, action='apply', display='',
                 freqcutoff=3.0, timecutoff=2.0,
                 combinescans=True, intent='*TARGET*',
                 winsize=1, flagdimension='time')
        flagdata(vis=msname, scan=itss, mode='tfcrop',
                 datacolumn='corrected', spw=cont_spw,
                 extendflags=False, action='apply', display='',
                 freqcutoff=3.0, timecutoff=2.0,
                 combinescans=True, intent='*TARGET*',
                 winsize=1, flagdimension='freq')
    # grow only
    for itss in idc_target_scans_str:
        flagdata(vis=msname, scan=itss, mode='extend',
                 datacolumn='corrected', spw=spw,
                 extendflags=True, action='apply', display='',
                 freqcutoff=100.0, timecutoff=100.0,
                 growtime=50.0, growfreq=30.0, intent='*TARGET*',
                 extendpols=False, combinescans=False,
                 winsize=1, flagdimension='freq')


def everythings_statwt(msname='', cont_spw='8:0~1024;3072~4096'):
    everythings_log("statwt")
    # Calibrators
    statwt(vis=msname, dorms=False, minsamp=2, intent='*CALIBRATE*',
           datacolumn='corrected')
    # Target
    statwt(vis=msname, dorms=False, fitspw=cont_spw, minsamp=2,
           intent='*TARGET*', datacolumn='corrected')
