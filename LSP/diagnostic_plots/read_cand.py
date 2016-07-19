import sys
import pandas as pd

OBS = sys.argv[1]
CAND = int(sys.argv[2])

WKDIR = '/projects/0/lotaas2/data/out/check/'+OBS

cands = pd.read_hdf(WKDIR+'/SinglePulses.hdf5','candidates')
cands = cands.loc[CAND]

print cands.DM

if cands.BEAM == 12: beamtype = 'incoherentstokes'
else: beamtype = 'stokes'

print beamtype

