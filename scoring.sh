#!/bin/bash

if [ $# -ne 1 ] ;then
    echo "This script creates a file with all the positive candidates within a pointing.It writes it in two locations:one locally inside the observation directory at the same level as the SAPs, and one on cartesius for the main log."
    echo ""
    echo "Usage:.sh <observation>"
    echo "i.e. arff_merger.sh L1212121"
    exit
fi

#SAPs
dirlist="SAP0 SAP1 SAP2"

#input data directory
obsdir="/projects/0/lotaas/data/out"
currentdir=`pwd`
cd $obsdir
absobsdir=`pwd -P`
cd $currentdir

#pointing name
obs="$1"

#predictions directories
predirlocal="$obsdir/$obs"
predircartesius="/projects/0/lotaas/data/predictions"
cd $predircartesius
abspredircartesius=`pwd -P`
cd $currentdir

#creating the predictions files
predictlocal="$predirlocal/predictions_$obs.txt"
predictcartesius="$predircartesius/predictions_$obs.txt"
rm $predictlocal $predictcartesius

touch $predictlocal
touch $predictcartesius
if [ ! -f $predictlocal ] ;then
    echo "$predictlocal file cannot be created."
    exit
fi
if [ ! -f $predictcartesius ] ;then
    echo "$predictcartesius file cannot be created."
    exit
fi

#get the absolute path
abspath=`pwd -P`

for sap in $dirlist; do
    #cd $obsdir/$obs/$sap
    for beam in `seq 0 73`; do   #was 72
	if [ -e $obsdir/$obs/$sap/$obs"_"$sap"_BEAM"$beam ];then
	    #cd $obs"_"$sap"_BEAM"$beam
	    if [ -f $obsdir/$obs/$sap/$obs"_"$sap"_BEAM"$beam/dmprof_$beam.arff ];then 
		    java -jar $abspath/PulsarProcessingScripts-master/ML.jar -v -m$abspath/classifierModels/DT_LOTAAS.model -o$abspredircartesius/temp-$1-$beam.txt -p$absobsdir/$obs/$sap/$obs"_"$sap"_BEAM"$beam/dmprof_$beam.arff -a1 >classifier.check
		    if [[ `grep -F "Could not make predictions" classifier.check` ]]
		    then
			echo "Classifier failed in SAP$sap""_""BEAM$beam">>$obsdir/$obs/classifier_fail_info
		    fi
		    rm classifier.check
		    if [ -f $abspredircartesius/temp-$1-$beam.txt ];then
			cat $abspredircartesius/temp-$1-$beam.txt >>$predictlocal
			cat $abspredircartesius/temp-$1-$beam.txt >>$predictcartesius
			rm $abspredircartesius/temp-$1-$beam.txt #cleaning up
		    else
			echo "The classifier did not create any output for SAP$sap""_""BEAM$beam" >>$obsdir/$obs/classifier_fail_info
		    fi
            else
		echo "$obsdir/$obs/$sap/$obs""_$sap""_BEAM$beam/dmprof_$beam.arff does not exist."
		echo "$obsdir/$obs/$sap/$obs""_$sap""_BEAM$beam/dmprof_$beam.arff" >>$predirlocal/absent.txt
		echo "$obsdir/$obs/$sap/$obs""_$sap""_BEAM$beam/dmprof_$beam.arff" >>$predircartesius/absent.txt
	    fi
	    #cd ..
	else
	    echo "$obsdir/$obs/$sap/$obs""_$sap""_BEAM$beam does not exist."
	    echo "$obsdir/$obs/$sap/$obs""_$sap""_BEAM$beam" >>$predirlocal/absent.txt
	    echo "$obsdir/$obs/$sap/$obs""_$sap""_BEAM$beam" >>$predircartesius/absent.txt
	fi
    done
done



