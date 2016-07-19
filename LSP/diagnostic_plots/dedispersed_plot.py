'''
Parse script to comunicate between dedispersed_plot.sh and LSPplot.timeseries
'''

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Parameters import *

OBS = sys.argv[1]
SAP = int(sys.argv[2])
BEAM = int(sys.argv[3])
CAND = int(sys.argv[4])
#INDIR = sys.argv[5]

WKDIR = '/projects/0/lotaas2/data/out/check/'+OBS


#mpl.rc('font',size=5)

cands = pd.read_hdf(WKDIR+'/SinglePulses.hdf5','candidates')
cands = cands.loc[CAND]

pulses = pd.read_hdf(WKDIR+'/SinglePulses.hdf5','pulses')
pulses = pulses[pulses.Candidate == CAND]
pulses.sort('Time',inplace=True)

nPlotBins = 21
nProfBins = 5

m = 1. / (4148.808/RES/4. * (F_MIN**-2 - F_MAX**-2) )  #Expected pulse inclination

for i,puls in enumerate(pulses.iterrows()):
  plt.clf()
  #plt.figure(figsize=(6,4))

  bin_peak = int(puls[1]['Sample'])
  DM_peak = puls[1]['DM']
  if DM_peak < 40.52: 
    duration = int(np.round(puls[1]['Duration']/0.00049152))
    DM_res = 0.01
  elif DM_peak < 141.77: 
    duration = int(np.round(puls[1]['Duration']/0.00049152/2.))
    DM_res = 0.05
  else: 
    duration = int(np.round(puls[1]['Duration']/0.00049152/4.))
    DM_res = 0.1

  nDMs = int(np.round(puls[1]['dDM'] / DM_res * 3))
  nBins = nPlotBins * nProfBins
  data = np.zeros((nDMs,nBins))

  scrunch_fact = int(np.round(duration / float(nProfBins)))
  if scrunch_fact < 1: scrunch_fact = 1

  bin_start = bin_peak - nBins/2 * scrunch_fact

  x = np.arange(nBins)
  DM_range = DM_peak - (nDMs/2 - np.arange(nDMs)) * DM_res 
  for j,DM in enumerate(DM_range):
    try:
      filename = '/timeseries/{0}_SAP{1}_BEAM{2}_DM{3:.2f}.dat'.format(OBS, SAP, BEAM, DM)
      ts = np.memmap(WKDIR + filename, dtype=np.float32, mode='r', offset=bin_start*4, shape=(nBins*scrunch_fact,))
      ts = np.mean(np.reshape(ts, (nBins, scrunch_fact)), axis=1)
      #ts -= np.median(ts)     #Set all the DM series at
      #ts /= np.abs(ts).max()  #the same flux range
      #ts *= DM_res * 2.  #Scale the distance between DM series
    except IOError: ts = np.zeros(nBins) + np.nan

    data[j] = ts
 
  #Line plot
  #for ts in data:
  #  plt.plot(x, ts+DM, 'k')
  #plt.plot(x, data[data.size/2+1], 'k', lw=2.)  #plot thick line for peak DM

  #Image plot
  plt.imshow(data,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=[0,nBins,DM_range[0], DM_range[-1]])

  #plt.ylim([DM_range[0]-DM_res, DM_range[-1]+DM_res*2])
  #plt.xlim([0,nBins])
  plt.xticks(x[::nBins/6], (x[::nBins/6]+bin_start)%1000, )
  plt.xlabel('Sample - {}000'.format(bin_start/1000))
  plt.ylabel('DM (pc/cc)')
  plt.title('{obs} SAP{sap} BEAM{beam} - Candidate {cand} Pulse {puls}'.format(obs=OBS,sap=SAP,beam=BEAM,cand=CAND,puls=i))

  #Plot contours  #FINIRE!!
  #y = -m*x
  #plt.plot(x,y+DM_range[-1],'r')

  #Inset profile
  filename = '/timeseries/{0}_SAP{1}_BEAM{2}_DM{3:.2f}.dat'.format(OBS, SAP, BEAM, DM_peak)
  ts = np.memmap(WKDIR + filename, dtype=np.float32, mode='r', offset=bin_start*4, shape=(nBins*scrunch_fact,))
  ts = np.mean(np.reshape(ts, (nBins, scrunch_fact)), axis=1)
  ts -= np.median(ts)
  #ts /= np.abs(ts).max()
  #ts *= DM_res * 2.
  
  plt.axes([.7, .7, .15, .15])
  plt.plot(x[nBins/2-scrunch_fact:nBins/2+scrunch_fact], ts[nBins/2-scrunch_fact:nBins/2+scrunch_fact]+DM_peak, 'k')
  plt.xticks([])
  plt.yticks([])

  plot_name = WKDIR+'/diagnostic_plots/{cand}/{cand}_{puls}_dedispersion.png'.format(cand=cands.id,puls=i)
  plt.savefig(plot_name,format='png',bbox_inches='tight',dpi=200)
 

