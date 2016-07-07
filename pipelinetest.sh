#!/bin/bash

#available modes
#copy: copies the raw data from the discoveries directory to raw
#delraw: removes the raw data from the raw directory
#delres: removes the processed results
#jobs: creates the jobs script
#process: runs the classifier and creates a txt file with the positives at the same level as the arff file

MODE=$1

#discoveries directory
DSCVR=/projects/0/lotaas/data/raw/discoveries
#raw data directory
RAW=/projects/0/lotaas/data/raw
#processed results directory
OUT=/projects/0/lotaas/data/out

#pick up all the discoveries (except a few)
PSRLIST=`ls ${DSCVR}`
PSRLISTCLEAN=`echo "${PSRLIST}"|grep -v "J1342+65_conf\|J1849+15_RRAT\|avg_profs"`

#pick up the obsid directories
OBSIDS=""
OBSIDDIRS=""
for i in ${PSRLISTCLEAN};do
    obsid=`ls ${DSCVR}/${i}`
    OBSIDS+="${obsid} "
    OBSIDDIRS+="${i}/${obsid} "
done

#copy the directories for processing
if [[ ${MODE} == "copy" ]];then
for i in ${OBSIDDIRS};do
    cp -r ${DSCVR}/${i} ${RAW}
done
exit
fi



#creating the batch processing script

if [[ ${MODE} == "jobs" ]];then

if [ -f ~/jobscript.sh ];then
    rm ~/jobscript.sh
fi

echo "#!/bin/bash" >~/jobscript.sh
postfix="_red"

for i in ${OBSIDS};do
    if [ -d "${RAW}/${i}/stokes" ];then
	STATE="stokes"
    elif [ -d "${RAW}/${i}/incoherentstokes" ];then
	STATE="incoherentstokes"
    elif [ -d "${RAW}/${i}/stokes" && -d "${RAW}/${i}/incoherentstokes" ];then
	STATE="stokes incoherentstokes"
    fi
    for s in ${STATE};do
	SAPS=`ls ${RAW}/${i}/${s}`
	for j in ${SAPS};do 
	    BEAMS=`ls ${RAW}/${i}/${s}/${j}`
	    for k in ${BEAMS};do
		echo "sbatch /home/sanidas/LOTAAS-Scripts/lotaas_fil_LTA_dips1.sh ${i%$postfix} ${j:3:1} ${k:4:2}" >>~/jobscript.sh
		echo "sleep 10" >>~/jobscript.sh
	    done
	done
    done
done
exit
fi

#post processing

if [[ ${MODE} == "process" ]];then

if [ -f ~/postprocessing.sh ];then
    rm ~/postprocessing.sh
fi

postfix="_red"
prefix="BEAM"
echo "#!/bin/bash" >~/postprocessing.sh

for i in ${OBSIDS};do
        SAPS=`ls ${OUT}/${i%$postfix}`
        for j in ${SAPS};do
            BEAMS=`ls ${OUT}/${i%$postfix}/${j}|grep "${i%$postfix}_
${j}_BEAM"`
            for k in ${BEAMS};do
		beamno=`echo "${k}"|awk -F_ '{print $3}'`
                echo "java -jar /home/sanidas/scripts/scores_robLyon/PulsarProcessingScripts-master/ML.jar -v -m/nfs/home6/sanidas/classifierModels/DT_LOTAAS.model -o/lustre1/0/lotaas/data/out/${i%$postfix}/${j}/${k}/classifier_results.txt -p/lustre1/0/lotaas/data/out/${i%$postfix}/${j}/${k}/dmprof_${beamno#$prefix}.arff -a1;sed -i 's/36scrunch/ps/g' /lustre1/0/lotaas/data/out/${i%$postfix}/${j}/${k}/classifier_results.txt" >>~/postprocessing.sh
            done
        done
    done

fi


#remove raw data directories

if [[ ${MODE} == "delraw" ]];then
    for i in ${OBSIDS};do
	rm -rf ${RAW}/${i}
    done
exit
fi

#remove results directories

postfix="_red"

if [[ ${MODE} == "delres" ]];then
    for i in ${OBSIDS};do
	rm -rf ${OUT}/${i%$postfix}
    done  
exit
fi


