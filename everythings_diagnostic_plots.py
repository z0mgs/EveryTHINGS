# General imports
import os
import numpy as np
# CASA imports
from tasks import plotms, plotcal, flagdata, flagcmd
from casac import casac

# CASA tools
msmd = casac.msmetadata()


def everythings_diagnostic_plots(msname='', plot_spw=8):
    # Find restfreq
    msmd.open(msname)
    freqs = msmd.chanfreqs(plot_spw)
    restfreq = str(np.mean(freqs) / 1000000) + 'MHz'
    del freqs
    msmd.close()
    plot_spw = str(plot_spw)
    #
    spath = 'EveryTHINGS_Diagnostic_Plots/'
    if not os.path.isdir(spath):
        os.mkdir(spath)
    #
    # amp-time
    plotms(vis=msname,
           xaxis='time', yaxis='amp', avgchannel='4096', coloraxis='ant1',
           avgantenna=True, overwrite=True, ydatacolumn='corrected',
           plotfile=spath+'amp-time_target.png', avgtime='10',
           showgui=False, intent='*TARGET*', spw=plot_spw)
    #
    plotms(vis='finalcalibrators.ms',
           xaxis='time', yaxis='amp', avgchannel='4096', coloraxis='ant1',
           avgantenna=True, overwrite=True, ydatacolumn='corrected',
           plotfile=spath+'amp-time_calibrators.png', avgtime='10',
           showgui=False, spw=plot_spw)
    #
    # amp-vel
    plotms(vis=msname,
           xaxis='vel', yaxis='amp', avgchannel='', coloraxis='ant1',
           avgtime='100000', avgscan=True, restfreq=restfreq,
           veldef='RADIO',
           avgantenna=True, overwrite=True, ydatacolumn='corrected',
           plotfile=spath+'amp-vel_target.png', intent='*TARGET*',
           showgui=False, exprange='all', spw=plot_spw)
    plotms(vis='finalcalibrators.ms',
           xaxis='vel', yaxis='amp', avgchannel='', coloraxis='ant1',
           avgtime='100000', avgscan=True, restfreq=restfreq,
           veldef='RADIO',
           avgantenna=True, overwrite=True, ydatacolumn='corrected',
           plotfile=spath+'amp-vel_calibrators.png',
           showgui=False, iteraxis='field', exprange='all', spw=plot_spw)
    #
    # phase-vel
    plotms(vis='finalcalibrators.ms',
           xaxis='vel', yaxis='phase', coloraxis='ant1',
           avgtime='100000', avgscan=True, restfreq=restfreq,
           veldef='RADIO',
           avgantenna=True, overwrite=True, ydatacolumn='corrected',
           plotfile=spath+'phase-vel.png',
           showgui=False, iteraxis='field', exprange='all',
           plotrange=[0, 0, -180, 180], spw=plot_spw)
    #
    # uvdist-amp
    plotms(vis=msname,
           xaxis='uvdist', yaxis='amp', avgchannel='4096', coloraxis='ant1',
           avgtime='5000', avgscan=False,
           avgantenna=False, overwrite=True, ydatacolumn='corrected',
           plotfile=spath+'amp-uv_target.png',
           showgui=False, iteraxis='field', exprange='all', intent='*TARGET*',
           spw=plot_spw)
    #
    plotms(vis='finalcalibrators.ms',
           xaxis='uvdist', yaxis='amp', avgchannel='4096', coloraxis='field',
           avgtime='5000', avgscan=False,
           avgantenna=False, overwrite=True, ydatacolumn='corrected',
           plotfile=spath+'amp-uv_calibrators.png',
           showgui=False, spw=plot_spw)
    #
    # uvdist-phase
    plotms(vis='finalcalibrators.ms',
           xaxis='uvdist', yaxis='phase', avgchannel='4096', coloraxis='ant1',
           avgtime='5000', avgscan=False,
           avgantenna=False, overwrite=True, ydatacolumn='corrected',
           plotfile=spath+'phase-uv.png',
           showgui=False, iteraxis='field', exprange='all',
           plotrange=[0, 0, -180, 180], spw=plot_spw)
    #
    # BPcal
    tblname = msname + '.hifv_finalcals.s13_4.finalBPcal.tbl'
    plotcal(caltable=tblname, xaxis='freq',
            yaxis='amp', iteration='antenna',
            subplot=651, overplot=False, clearpanel='Auto', showflags=False,
            showgui=False, figfile=spath+'finalBPcal.png', spw=plot_spw)
    #
    # UV coverage
    plotms(vis=msname, intent='*TARGET*', xaxis='Uwave', yaxis='Vwave',
           avgchannel='4096', coloraxis='ant1', customflaggedsymbol=True,
           showgui=False, plotfile=spath+'uvcoverage.png', spw=plot_spw)
    #
    # flagdata summary
    flagsummary = flagdata(msname, mode='summary', spw=plot_spw)
    antsummary = flagsummary['antenna']
    flaglist = []
    for ant in antsummary.keys():
        flaglist.append([ant, round(antsummary[ant]['flagged'] /
                                    antsummary[ant]['total'], 3)])
    flaglist = np.sort(flaglist, axis=0)
    np.savetxt(spath + 'flagsummary_all.txt', flaglist, fmt='%s')
    #
    flagsummary = flagdata(msname, mode='summary', intent='*CALI*',
                           spw=plot_spw)
    antsummary = flagsummary['antenna']
    flaglist = []
    for ant in antsummary.keys():
        flaglist.append([ant, round(antsummary[ant]['flagged'] /
                                    antsummary[ant]['total'], 3)])
    flaglist = np.sort(flaglist, axis=0)
    np.savetxt(spath+'flagsummary_cali.txt', flaglist, fmt='%s')
    #
    flagcmd(msname, inpmode='list', inpfile=msname[:-3] + '.flagcmds.txt',
            action='plot', plotfile=spath + 'onlineFlags.png')
