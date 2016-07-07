#!/bin/bash

if [ $# -ne 1 ] ;then
    echo "This script creates a file with all the positive candidates within a pointing.It writes it in two locations:one locally inside the observation directory at the same level as the SAPs, and one on cartesius for the main log."
    echo ""
    echo "Usage:.sh <observation>"
    echo "i.e. arff_merger.sh L1212121"
    exit
fi

#SAPs
dirlist="SAP0"

#input data directory
obsdir="/projects/0/lotaas2/data/out/confirmations"

#pointing name
obs="$1"

#predictions directories
predirlocal="$obsdir/$obs"
predircartesius="/projects/0/lotaas2/data/predictions"

#the merged arff file
#arffmerged="$obsdir/$1/dmprof_$1""_merged.arff"
#touch $arffmerged
#if [ ! -f $arffmerged ] ;then
#    echo "Merged arff file cannot be created."
#    exit
#fi
#rm $arffmerged #security measure in case this script is executed manually more than once 

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

currdir=`pwd`
cd $LOTAAS_PL
truedir=`pwd -P`
cd ${predircartesius}
truedirpred=`pwd -P`
cd ${obsdir}
trueobsdir=`pwd -P`
cd ${currdir}

for sap in $dirlist; do
    #cd $obsdir/$obs/$sap
    for beam in `seq 0 127`; do   
	if [ -e $obsdir/$obs/$sap/$obs"_"$sap"_BEAM"$beam ];then
	    #cd $obs"_"$sap"_BEAM"$beam
	    if [ -f $obsdir/$obs/$sap/$obs"_"$sap"_BEAM"$beam/dmprof_$beam.arff ];then 
		    java -jar $LOTAAS_PL/scores_robLyon/PulsarProcessingScripts-master/ML.jar -v -m${truedir}/classifierModels/DT_LOTAAS.model -o${truedirpred}/temp-$1-$beam.txt -p${trueobsdir}/$obs/$sap/$obs"_"$sap"_BEAM"$beam/dmprof_$beam.arff -a1
		    if [ -f ${truedirpred}/temp-$1-$beam.txt ];then
			sed -i 's|./|/projects/0/lotaas2/data/out/confirmations/'$obs/$sap/$obs"_"$sap"_BEAM"$beam/"BEAM"$beam"_sift/"'|' ${predircartesius}/temp-$1-$beam.txt
			cat ${truedirpred}/temp-$1-$beam.txt >>$predictlocal
			cat ${truedirpred}/temp-$1-$beam.txt >>$predictcartesius
			rm ${truedirpred}/temp-$1-$beam.txt #cleaning up
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



