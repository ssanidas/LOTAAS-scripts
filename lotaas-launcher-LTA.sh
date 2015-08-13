#!/bin/bash

# launch script for multiple Cartesius processing (all beams + saps in pointing)

# $1 is pointing e.g, L155450

doProcess=false

for sap in {0..2}
do
    for beam in {0..73}
    do	
	
	# check if processed beam dir exists and contains pfd files.
	if [ -e /projects/0/lotaas/data/out/$1/SAP$sap/$1_SAP${sap}_BEAM${beam}/BEAM${beam}_sift/FOLDS/CORE_0 ]
	then
	if [ `ls /projects/0/lotaas/data/out/$1/SAP$sap/$1_SAP${sap}_BEAM${beam}/BEAM${beam}_sift/FOLDS/CORE_0/*.pfd* | wc -l` -gt 1 ]
	then		
	   echo "pfd's exist"
	   continue
	fi
        fi
	
	if [ ${beam} -eq 12 ]
	then
	    stokes="incoherentstokes"
	else
	    stokes="stokes"
	fi

        if [ -e /projects/0/lotaas/data/raw/$1_red/$stokes/SAP${sap}/BEAM${beam}/*_SAP${sap}_BEAM${beam}.fits ] 
	then
            if [ -s /projects/0/lotaas/data/raw/$1_red/$stokes/SAP${sap}/BEAM${beam}/*_SAP${sap}_BEAM${beam}.fits ]
	    then
		
		while [ `squeue -u sanidas | grep lotaas | wc -l` -gt 150 ]
		do
		    echo 'sleeping for 15 seconds'
		    sleep 20
		done
		
			
		sbatch /home/sanidas/LOTAAS-Scripts/lotaas_fil_LTA.sh $1 ${sap} ${beam} >job_id_$1.txt
		cat job_id_$1.txt|awk '{print $4}' >>jobs_submitted_$1.txt
		rm job_id_$1.txt
		sleep 20
		
	    else
		echo "Fits file $1_SAP${sap}_BEAM${beam}.fits is zero size"
	    fi
	else
	    echo "$1_SAP${sap}_BEAM${beam}.fits is not present"
	fi
        
	
    done
done

squeue -u sanidas >jobsonqueue.txt
jobcheck=`grep -f jobs_submitted_$1.txt jobsonqueue.txt|wc -l`
while [[ $jobcheck != 0 ]];do
    echo "Waiting for all jobs to finish"
    sleep 600
    squeue -u sanidas >jobsonqueue.txt
    jobcheck=`grep -f jobs_submitted_$1.txt jobsonqueue.txt|wc -l`
done

#SP analysis
#cd /home/sanidas/scripts_dev
sh LSPs_launcher.sh $1

#candidates
#cd /home/sanidas/scripts_dev
sh scoring.sh $1
sh predictions.sh $1

#cleaning up
cd /home/sanidas/LOTAAS-Scripts
rm jobs_submitted_$1.txt

