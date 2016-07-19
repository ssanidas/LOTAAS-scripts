'''
Parse script to comunicate between store_online.sh and Internet.py
'''

from Internet import upload
import sys
import pandas as pd

OBS = sys.argv[1]

WKDIR = '/projects/0/lotaas2/data/out/check/'+OBS

cands = pd.read_hdf(WKDIR+'/SinglePulses.hdf5','candidates')
cands = cands[cands.main_cand==0].head(30)

folder = WKDIR+'/products/candidates/.'

#try: 
upload(cands,OBS,folder)
#  print "Observation stored"
#except: 
#  print "ATTENTION!\n\nConnession problem, update candidates in a second moment\n\n"
#  with open('/home/danielem/LOTAAS_ERRORS.txt','a') as f:
#    f.write("ATTENTION! Connession problemwith obs. {} while restoring\n".format(OBS))



