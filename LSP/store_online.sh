#!/bin/bash

# Store online a processed pointing again.
# USAGE: sh store_online.sh OBS_ID

if [ $# -ne 1 ] ;then
    echo USAGE: store_online.sh OBS_ID
    exit
fi

OBS=$1
WRKDIR=/projects/0/lotaas2/data/out/check/${OBS}
mkdir $WRKDIR



if [ ! -f ${WRKDIR}/SinglePulses.hdf5 ]; then
  #fetch the tarball
  cd /archive/pulsar/lotaas/out_tarballs
  echo "Tarballs:"
  echo `dmls -l ${OBS}_sp* | awk '{print $9}' | tail -1`
  arcfile=`dmls -l ${OBS}_sp* | awk '{print $9}' | tail -1`
  dmget -a $arcfile
  tar -xf $arcfile -C ${WRKDIR} sp/products/ --strip-components 1 
fi

cd $WRKDIR

python /projects/0/lotaas2/software/LSP/store_online.py ${OBS}

