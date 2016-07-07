#!/bin/bash

# launch script for multiple Cartesius processing (all beams + saps in pointing)

# $1 is pointing e.g, L155450


rawdir=/projects/0/lotaas2/data/raw/confirmations
outdir=/projects/0/lotaas2/data/out/confirmations

sap=0
for beam in {0..127}
do	
	
    if [ ${beam} -eq 0 ]
    then
	stokes="incoherentstokes"
    else
	stokes="stokes"
    fi
    
    if [ -e ${rawdir}/$1_red/$stokes/SAP${sap}/BEAM${beam}/*_SAP${sap}_BEAM${beam}.fits ] 
	then
        if [ -s ${rawdir}/$1_red/$stokes/SAP${sap}/BEAM${beam}/*_SAP${sap}_BEAM${beam}.fits ]
	then
	    
	    while [ `squeue -u ${USER} | grep lotaas | wc -l` -gt 150 ]
	    do
		echo 'sleeping for 30 seconds'
		sleep 15
	    done
		
	    
	    sbatch $LOTAAS_PL/lotaas_fil_confirmations_nosp_msp.sh $1 ${sap} ${beam}  >job_id_$1.txt
                cat job_id_$1.txt|awk '{print $4}' >>jobs_submitted_$1.txt
                rm job_id_$1.txt

	    sleep 30
			
	else
	    echo "Fits file $1_SAP${sap}_BEAM${beam}.fits is zero size"
	fi
    else
	echo "$1_SAP${sap}_BEAM${beam}.fits is not present"
    fi
    
    
done

squeue -u ${USER} >jobsonqueue.txt
jobcheck=`grep -f jobs_submitted_$1.txt jobsonqueue.txt|wc -l`
while [[ $jobcheck != 0 ]];do
    echo "Waiting for all jobs to finish"
    sleep 600
    squeue -u ${USER} >jobsonqueue.txt
    jobcheck=`grep -f jobs_submitted_$1.txt jobsonqueue.txt|wc -l`
done


for i in `seq 0 127`;do
    cp ${outdir}/$1/SAP0/dmprof_$i.arff ${outdir}/$1/SAP0/$1"_"SAP0_BEAM$i
done


sh $LOTAAS_PL/run_classifier_confirmations.sh $1
sh $LOTAAS_PL/predictions_confirmations.sh $1

#heatmap plot
cd ${outdir}/$1
parfile=`ls ${rawdir}/${1}_red/*.par`
candname=`cat ${parfile}|grep PSR|awk "{print $2}"`
candDM=`cat ${parfile}|grep dm|awk "{print $2}"|awk -F. "{print $1}"`
sh $LOTAAS_PL/heatmap_conf.sh ${candDM} ${candname} $1

#cleaning up                                                                                               

rm jobs_submitted_$1.txt
