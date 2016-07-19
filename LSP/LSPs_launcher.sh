#!/bin/bash -x


#####################################
# 
# LOTAAS Single Pulse search launcher
#
# written by Daniele Michilli
# 
# usage: LSPs_launcher <OBS_ID>
#
#####################################


#OBS_PATH=/projects/0/lotaas/data/out/new/

OBS_PATH=/projects/0/lotaas2/data/out/
#OBS_PATH=/projects/0/lotaas/data/out/new/archived/
#OBS_PATH=/projects/0/lotaas/data/out/confirmations/new/

ind=0

for i; do

  if [ -e $OBS_PATH/$i ]; then
    if [ ! -e $OBS_PATH/$i/sp ]; then
      if [ `ls $OBS_PATH/$i/SAP*/$i*/BEAM*_sift/sp/*_singlepulse.tgz | wc -l` -gt 0 ]; then
        echo "Job submitted for obs. $i"
        echo "Job submitted for obs. $i" #>> obs_submitted.log
        sbatch $LOTAAS_PL/../LSP/LSPs_job.sh $i #$OBS_PATH
      else
        echo "Obs. $i doesn't contain any Single Pulse file"
      fi
    else
      echo "Obs. $i has already been analyzed"
      echo "Remove the folder $OBS_PATH/$i/sp to perform a new analysis"
    fi
  else
    echo "Obs. $OBS_PATH/$i doesn't exist"
  fi

done

exit 0

