#!/bin/bash


#need to run from same direcrtory as ML.jar

#LOCATION OF ARFF SCORES /projects/lotaas/data/out/$L_point/dmprof_$L_point.arff

Lfile=/home/sanidas/scripts/new_processed.lis

for L_point in $(cat $Lfile)
do
    if [ -s /projects/0/lotaas/data/out/$L_point/dmprof_13.arff ]
    then
	echo 'predicting dmprof_$L_point.arff'
	java -jar ML.jar -v -m/nfs/home6/sanidas/classifierModels/DT_LOTAAS.model -o/lustre1/0/lotaas/data/periodicity/predict.txt -p/lustre1/0/lotaas/data/out/$L_point/dmprof_13.arff -a1
    else
	echo '$L_point/dmprof_$L_point.arff does not exist'
    fi
done

