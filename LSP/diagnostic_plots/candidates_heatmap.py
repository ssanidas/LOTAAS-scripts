'''
Parse script to comunicate between candidates_heatmap.sh and LSPplot.heatmap 
'''

#from LSPplot import heatmap
#TEMP!
import imp
a = imp.load_source('b', '/home/danielem/LSPs/src/LSPplot.py')
from b import heatmap


import sys
import pandas as pd

OBS = sys.argv[1]
SAP = int(sys.argv[2])
BEAM = int(sys.argv[3])
CAND = int(sys.argv[4])

WKDIR = '/projects/0/lotaas2/data/out/check/'+OBS

cands = pd.read_hdf(WKDIR+'/SinglePulses.hdf5','candidates')
cands = cands.loc[CAND]

pulses = pd.read_hdf(WKDIR+'/SinglePulses.hdf5','pulses')
pulses = pulses[pulses.Candidate == CAND]
pulses.sort('Time',inplace=True)

events = pd.read_hdf(WKDIR+'/SinglePulses.hdf5','events')
events = events[(events.SAP == SAP) & (events.BEAM > 12)]
#events = events[(events.DM > cands.DM - 0.501) & (events.DM < cands.DM + 0.501)]

for i,n in enumerate(pulses.iterrows()):
  ev = events[(events.Time > n[1]['Time'] - n[1]['Duration']*2) & (events.Time < n[1]['Time'] + n[1]['Duration']*2) & (events.DM > n[1]['DM'] - 0.201) & (events.DM < n[1]['DM'] + 0.201)]
  plot_name = WKDIR+'/diagnostic_plots/{cand}/{obs}_{sap}_{beam}_{cand}_{puls}_heatmap.png'.format(cand=cands.id,obs=OBS,sap=SAP,beam=BEAM,puls=i)
  heatmap(ev,plot_name,idL=OBS,sap=SAP,beam=BEAM,cand=CAND,pulse=i,dm=n[1]['DM'],time=n[1]['Time'],duration=n[1]['Duration']) 
  
print "Plots stored in {}".format(WKDIR) 


