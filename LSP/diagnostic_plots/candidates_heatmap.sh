#!/bin/bash

# Produce a cumulative SNR heatmap for the beams within a SAP 
# around the pulses of a single-pulse candidate indicated by the user.
#
# USAGE: sh candidates_heatmap.sh candID, where candID is the first column of LSP_candidates spreadsheet
# EXAMPLE: sh candidates_heatmap.sh L101010_1_1_10

if [ $# -ne 1 ] ;then
    echo USAGE: sh candidates_heatmap.sh candID, where candID is the first column of LSP_candidates spreadsheet
    echo EXAMPLE: sh candidates_heatmap.sh L101010_1_1_10
    exit
fi

candID=$1

arr=(${candID//_/ })
OBS=${arr[0]}
SAP=${arr[1]}
BEAM=${arr[2]}
CAND=${arr[3]}

if [ $BEAM -le 12 ]; then
  echo "It is not possible to produce a heatmap for beam numbers lower than 13"
  exit
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
WRKDIR=/projects/0/lotaas2/data/out/check/${OBS}
mkdir $WRKDIR

if [ ! -f $WRKDIR/SinglePulses.hdf5 ]; then
  if [ -f /projects/0/lotaas2/data/out/${OBS}/sp/products/SinglePulses.hdf5 ]; then
    cp /projects/0/lotaas2/data/out/${OBS}/sp/products/SinglePulses.hdf5 $WRKDIR/SinglePulses.hdf5
  elif [ -f /projects/0/lotaas2/data/out/new/${OBS}/sp/products/SinglePulses.hdf5 ]; then
    cp /projects/0/lotaas2/data/out/new/${OBS}/sp/products/SinglePulses.hdf5 $WRKDIR/SinglePulses.hdf5
  else
    #fetch the tarball
    cd /archive/pulsar/lotaas/out_tarballs
    echo "Tarballs:"
    echo `dmls -l ${OBS}_sp* | awk '{print $9}' | tail -1`
    arcfile=`dmls -l ${OBS}_sp* | awk '{print $9}' | tail -1`
    dmget -a $arcfile
    tar -xf $arcfile -C ${WRKDIR} sp/products/SinglePulses.hdf5 --strip-components 2
  fi
fi

#Call the python script to produce the plots
cd $WRKDIR
mkdir diagnostic_plots
mkdir diagnostic_plots/${candID}
python ${SCRIPT_DIR}/candidates_heatmap.py ${OBS} ${SAP} ${BEAM} ${CAND}

scp -prq ${WRKDIR}/diagnostic_plots/${candID} ag004:/var/www/lofarpwg/lotaas-sp/observations/${OBS}

