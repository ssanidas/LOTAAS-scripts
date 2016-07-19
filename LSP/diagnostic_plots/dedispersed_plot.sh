#!/bin/bash

# Plot a range of 10 dedispersed timeseries centered around the candidate pulses
#
# USAGE: sh 
# EXAMPLE: sh 

if [ $# -ne 3 ] ;then
    echo USAGE: sh candidates_heatmap.sh candID, where candID is the first column of LSP_candidates spreadsheet
    echo EXAMPLE: sh candidates_heatmap.sh L101010_1_1_10
    exit
fi

candID=$1
#fileID=$2

arr=(${candID//_/ })
OBS=${arr[0]}
SAP=${arr[1]}
BEAM=${arr[2]}
CAND=${arr[3]}

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

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
arr=()
while read line ; do
  arr+=($line)
done < <(python ${SCRIPT_DIR}/read_cand.py ${OBS} ${CAND})
DM=${arr[0]}
beamtype=${arr[1]}

INDIR=/projects/0/lotaas/data/raw/${fileID}/${beamtype}/SAP${SAP}/BEAM${BEAM}

sh /projects/0/lotaas2/software/LSP/spdspsr.sh ${fileID} ${SAP} ${BEAM} ${DM} 1

python ${SCRIPT_DIR}/dedispersed_plot.py $OBS $SAP $BEAM $CAND $INDIR

scp -prq ${WRKDIR}/diagnostic_plots/${candID} ag004:/var/www/lofarpwg/lotaas-sp/observations/${OBS}

