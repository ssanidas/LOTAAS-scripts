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
obsdir="/projects/0/lotaas2/data/out"

#pointing name
obs="$1"

#predictions directories
predirlocal="$obsdir/$obs"
predircartesius="/projects/0/lotaas2/data/predictions"


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

currentdir=`pwd`
cd $LOTAAS_PL
realrepdir=`pwd -P`
cd ${predircartesius}
realpreddir=`pwd -P`
cd ${obsdir}
realobsdir=`pwd -P`
cd ${currentdir}

for sap in $dirlist; do
    #cd $obsdir/$obs/$sap
    for beam in `seq 0 73`; do  
	if [ -e $obsdir/$obs/$sap/$obs"_"$sap"_BEAM"$beam ];then
	    #cd $obs"_"$sap"_BEAM"$beam
	    if [ -f $obsdir/$obs/$sap/$obs"_"$sap"_BEAM"$beam/dmprof_$beam.arff ];then 
		    java -jar $LOTAAS_PL/scores_robLyon/PulsarProcessingScripts-master/ML.jar -v -m${realrepdir}/classifierModels/DT_LOTAAS.model -o${realpreddir}/temp-$1-$beam.txt -p${realobsdir}/$obs/$sap/$obs"_"$sap"_BEAM"$beam/dmprof_$beam.arff -a1 >classifier.check.${obs}
		    if [[ `grep -F "Could not make predictions" classifier.check.${obs}` ]]
		    then
			echo "Classifier failed in SAP${sap}_BEAM${beam}">>$obsdir/$obs/classifier_fail_info
		    fi
		    rm classifier.check.${obs}
		    if [ -f ${realpreddir}/temp-$1-$beam.txt ];then
			cat ${realpreddir}/temp-$1-$beam.txt >>$predictlocal
			cat ${realpreddir}/temp-$1-$beam.txt >>$predictcartesius
			rm ${realpreddir}/temp-$1-$beam.txt #cleaning up
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



