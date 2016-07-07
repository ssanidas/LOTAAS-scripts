from PIL import Image
import matplotlib.pyplot as plt
import matplotlib as mp
import matplotlib.image as mpimg
import numpy as Num
import parameters
import string
import time
from math import *
#from scipy.optimize import leastsq
import struct
#import sys, psr_utils
#import psr_utils
import sys
from types import StringType, FloatType, IntType
from bestprof import bestprof
import glob
import pylab
import shutil
import os


obs="L404256"

fig=plt.figure()

psr_good = open("/projects/0/lotaas/data/sotiris_viewing/goodNameOnly.lis",'a')
psr_ok = open("/projects/0/lotaas/data/sotiris_viewing/okNameOnly.lis",'a')
psr_known = open("/projects/0/lotaas/data/sotiris_viewing/knownNameOnly.lis", 'a')
rfi = open("/projects/0/lotaas/data/sotiris_viewing/rfiNameOnly.lis", "a")

list1 = open("/projects/0/lotaas/data/sotiris_viewing/knownNameOnly.lis", "r").readlines()
list2 = open("/projects/0/lotaas/data/sotiris_viewing/okNameOnly.lis", "r").readlines()
list3 = open("/projects/0/lotaas/data/sotiris_viewing/goodNameOnly.lis", "r").readlines()
list4 = open("/projects/0/lotaas/data/sotiris_viewing/rfiNameOnly.lis", "r").readlines()
list5 = open("/projects/0/lotaas/data/sotiris_viewing/alllongnameOnly.lis", "r").readlines()
list6 = open("/projects/0/lotaas/data/sotiris_viewing/knownsort.lis", "r").readlines()


list_all = list4 + list5 + list6 + list1 + list2 + list3 

good_dst = "/projects/0/lotaas/data/cands_sotiris/cands1/"
ok_dst = "/projects/0/lotaas/data/cands_sotiris/cands2/"
known_dst = "/projects/0/lotaas/data/cands_sotiris/known_pulsars/to_be_sorted/"

#for candpfd1 in open('/projects/lotaas/data/cands/view/candpfd.lis', "r").readlines(): 
#for candpfd1 in open('/projects/lotaas/data/cands/view/newrandom', "r").readlines(): 
for candpfd in glob.glob('/projects/0/lotaas/data/pcands/' + obs +'/*.png'):
    
    candpfd2 = candpfd.strip()
    candpfd1 = (candpfd2.split('/')[7]).split('.png')[0]
    print candpfd1
    #DM = (candpfd.split('DM')[1]).split('.')[0]
    
#    if DM == '26':
#        print candpfd + "skipping"
#        continue
#    if DM == '43':
#        print candpfd + "skipping"
#        continue
    
    if os.path.isfile(good_dst + candpfd1 + '.png'):
        print "in good dir" 
        continue 
    elif os.path.isfile(ok_dst + candpfd1 + '.png'):
        print "in ok dir"
        continue
    elif os.path.isfile(known_dst + candpfd1 + '.png'):
        print "known psr"
        continue
    skip = 0
    
    for old in list_all:
        if string.strip(old) == candpfd1:
            skip = 1
    
    if skip == 1:
        continue
    
        
    fig.set_tight_layout(True)
    plt.ion()
    image=mpimg.imread(candpfd)
    ax = fig.add_subplot(1, 1, 1)
    imgplot = plt.imshow(image)
    plt.show(block=False)
    
    
    a = raw_input('Next plot? 1(good), 2(ok), k(nown), r(fi), q(uit) \n')
    if a == "1":
        psr_good.write(candpfd1 + '\n')
        shutil.copy(candpfd, good_dst)
        
    if a == "2":
        psr_ok.write(candpfd1 + '\n')
        shutil.copy(candpfd, ok_dst)
        
    if a == "k":
        psr_known.write(candpfd1 + '\n')
        shutil.copy(candpfd, known_dst)
        
    if a == "r":
        rfi.write(candpfd1 + '\n')
    if a == "q":
        exit(-1)
    
    
    plt.clf()

plt.close()
