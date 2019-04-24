"""
Reference: Astroua's LocalGroup-VLA pipeline
We refer heavily to their method of customized bandpass calibration
https://github.com/Astroua/LocalGroup-VLA/blob/master/14A-235/
pipeline_scripts/casa_pipeline_lines.py
"""
import os
from copy import copy
from glob import glob
# CASA imports
import pipeline.infrastructure.casatools as casatools
import pipeline.hif.heuristics.findrefant as findrefant
from tasks import bandpass


# Post message to casa log and terminal, with EveryTHINGS marker
def everythings_log(msg):
    casatools.post_to_log("::EveryTHINGS:: " + msg)


# Make sure I get the correct table
def get_table_index(s):
    s2 = s.split('_')
    s3 = s2[1].split('.')
    s4 = s3[1].strip('s')
    return int(s4)


bpcal_tbl = {'test': "{}.hifv_testBPdcals.s*_4.testBPcal.tbl",
             'semifinal': "{}.hifv_semiFinalBPdcals.s*_4.BPcal.tbl",
             'final': "{}.hifv_finalcals.s*_4.finalBPcal.tbl"}

delay_tbl = {'test': "{}.hifv_testBPdcals.s*_2.testdelay.tbl",
             'semifinal': "{}.hifv_semiFinalBPdcals.s*_2.delay.tbl",
             'final': "{}.hifv_finalcals.s*_2.finaldelay.tbl"}

gain_tbl = {'test': "{}.hifv_testBPdcals.s*_3.testBPdinitialgain.tbl",
            'semifinal': "{}.hifv_semiFinalBPdcals.s*_3.BPinitialgain.tbl",
            'final': "{}.hifv_finalcals.s*_3.finalBPinitialgain.tbl"}

len_tbl = {'test': 1, 'semifinal': 2, 'final': 1}


def astroua_bandpass(myvis='', stage='test', context=None,
                     refantignore=''):
    everythings_log('Begin bandpass fillgap correction script.')
    everythings_log('Bandpass stage: ' + stage)
    # stage = 'test', 'semifinal', 'final'
    # Look for BP table
    bpname = glob(bpcal_tbl[stage].format(myvis))
    assert len(bpname) == len_tbl[stage]
    bpname.sort(reverse=True, key=get_table_index)
    bpname = bpname[0]
    assert os.path.isdir(bpname)
    everythings_log('bpname: ' + bpname)
    # Remove already-made version
    # rmtables(bpname)
    # Or copy to another name to check against
    os.system("mv {0} {0}.orig".format(bpname))
    # Get the scan/field selections
    scanheur = context.evla['msinfo'][context.evla['msinfo'].keys()[0]]
    everythings_log('field: ' + scanheur.bandpass_field_select_string)
    everythings_log('scan: ' + scanheur.bandpass_scan_select_string)
    # Grab list of preferred refants
    refantfield = scanheur.calibrator_field_select_string
    refantobj = findrefant.RefAntHeuristics(vis=myvis, field=refantfield,
                                            geometry=True, flagging=True,
                                            intent='',
                                            spw='',
                                            refantignore=refantignore)
    RefAntOutput = refantobj.calculate()
    refAnt = ','.join(RefAntOutput)
    everythings_log('refAnt: ' + refAnt)
    #
    # priorcals part
    #
    # Lastly get list of other cal tables to use in the solution
    gc_tbl = glob("{}.hifv_priorcals.s*_2.gc.tbl".format(myvis))
    assert len(gc_tbl) == 1
    opac_tbl = glob("{}.hifv_priorcals.s*_3.opac.tbl".format(myvis))
    assert len(opac_tbl) == 1
    rq_tbl = glob("{}.hifv_priorcals.s*_4.rq.tbl".format(myvis))
    assert len(rq_tbl) == 1
    swpow_tbl = glob("{}.hifv_priorcals.s*.swpow.tbl".format(myvis))
    # Check ant correction
    ant_tbl = glob("{}.hifv_priorcals.s*_6.ants.tbl".format(myvis))
    #
    priorcals = [gc_tbl[0], opac_tbl[0], rq_tbl[0]]
    everythings_log('gc_tbl: ' + gc_tbl[0])
    everythings_log('opac_tbl: ' + opac_tbl[0])
    everythings_log('rq_tbl: ' + rq_tbl[0])
    #
    if len(ant_tbl) == 1:
        priorcals.extend(ant_tbl)
    if len(swpow_tbl) == 1:
        priorcals.extend(swpow_tbl)
    everythings_log('ant_tbl: ' + ant_tbl[0])
    #
    # priorcals part ends
    #
    del_tbl = glob(delay_tbl[stage].format(myvis))
    assert len(del_tbl) == len_tbl[stage]
    del_tbl.sort(reverse=True, key=get_table_index)
    del_tbl = del_tbl[0]
    #
    BPinit_tbl = glob(gain_tbl[stage].format(myvis))
    assert len(BPinit_tbl) == len_tbl[stage]
    BPinit_tbl.sort(reverse=True, key=get_table_index)
    BPinit_tbl = BPinit_tbl[0]
    #
    gaintables = copy(priorcals)
    gaintables.extend([del_tbl, BPinit_tbl])
    everythings_log('del_tbl: ' + del_tbl)
    everythings_log('BPinit_tbl: ' + BPinit_tbl)
    for tbl in gaintables:
        assert os.path.isdir(tbl)
    #
    everythings_log('Begin bandpass calibration')
    bandpass(vis=myvis,
             caltable=bpname,
             field=scanheur.bandpass_field_select_string,
             selectdata=True,
             scan=scanheur.bandpass_scan_select_string,
             solint='inf',
             combine='scan',
             refant=refAnt,
             minblperant=4,
             minsnr=5.0,
             solnorm=False,
             bandtype='B',
             smodel=[],
             append=False,
             fillgaps=400,
             docallib=False,
             gaintable=gaintables,
             gainfield=[''],
             interp=[''],
             spwmap=[],
             parang=True)
