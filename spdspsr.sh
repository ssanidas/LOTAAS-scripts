#!/bin/bash

#================================#
#SP CANDIDATE DEDISPERSING SCRIPT#
#================================#

#Description:
#This script retrieves a particular beam from the raw data,
#and dedisperses the fits file to a specific DM and with a
#specific downsampling factor


if [ $# -ne 5 ];then
    echo ""
    echo "The scripts has to be given five arguments to execute: FileID,SAP,BEAM,DM,Downsampling Factor"
    echo ""
    echo "i.e., sh sp_cand_folding.sh L121212 2 37 4.87 2"
    exit
fi

#PRELIMINARIES
#raw data directory
RAWDIR="/projects/0/lotaas2/data/raw"
#output directory 
OUTDIR="/projects/0/lotaas2/data/out/check"

#Observation ID
OBSID=$1
OBSIDNUM=`echo ${OBSID}|sed "s/L//g"`
SAP=$2
BEAM=$3
DM=$4
DS=$5

if [ ! -d ${RAWDIR}/${OBSID} ];then
    
    #Create the csv file of the project
    echo "================================"
    cd ${RAWDIR}
    echo "Creating the appropriate csv file..." 
    if [[ "${OBSIDNUM}" -lt 35000 && "${OBSIDNUM}" -gt 25000 ]];then
	PRCODE="LC3_014"
    elif [[ "${OBSIDNUM}" -lt 403000 && "${OBSIDNUM}" -gt 344621 ]];then
	PRCODE="LC4_030"
    elif [[ "${OBSIDNUM}" -lt 249928 ]];then
	echo "The requested ObsID could not be located in any of the LC3_014, LC4_030 or LT5_004 data. Exiting..."
	exit
    else
	PRCODE="LT5_004"
    fi

    echo "Project code:${PRCODE}"
    echo "Creating the csv file of ${PRCODE}..."
#    ./home/sanidas/scripts_lotaas/lta-query.py -p ${PRCODE}
    echo ".csv file of ${PRCODE} created."
    
    #Download the particular beam
    CSVFN=`echo "${PRCODE}.csv" | awk '{print tolower($0)}'` #get the low-case version of the csv file name

    $LOTAAS_PL/lta-retrieve.py --csvfile=${CSVFN} --sap=${SAP} --tab=${BEAM} ${OBSID}
fi

#dedispersing - ATTENTION! it will need modification for early pointings
if [[ ${BEAM} -eq "12" ]];then
    beamtype="incoherentstokes"
else
    beamtype="stokes"
fi
    
cd ${RAWDIR}/${OBSID}/${beamtype}/SAP${SAP}/BEAM${BEAM}
FITSFNAME=`ls *.fits|awk -F. '{print $1}'`
STAMP=`date|awk '{print $4}'|sed 's/\://g'`

if [ ! -d ${OUTDIR}/${OBSID} ];then
    mkdir -p ${OUTDIR}/${OBSID}
fi

#Creating the fil file

if [ ! -f ${RAWDIR}/${OBSID}/${beamtype}/SAP${SAP}/BEAM${BEAM}/${FITSFNAME}.fil ];then
    python $LOTAAS_PL/psrfits2fil.py --noscales --noweights --nooffsets ${FITSFNAME}.fits ${FITSFNAME}.fil
fi

#Recreating the RFI mask
if [ ! -f ${RAWDIR}/${OBSID}/${beamtype}/SAP${SAP}/BEAM${BEAM}/${OBSID}_SAP${SAP}_BEAM${BEAM}_spdspsr_rfifind.mask ];then
    rfifind -blocks 10 -noweights -noscales -nooffsets -o ${RAWDIR}/${OBSID}/${beamtype}/SAP${SAP}/BEAM${BEAM}/${OBSID}_SAP${SAP}_BEAM${BEAM}_spdspsr *.fil
fi

#prepfold -nsub 288 -noxwin -ndmfact 1 -p ${PERIOD} -dm ${DM} -o ${OUTDIR}/${OBSID}/manual_fold_${STAMP} -n 100 -pstep 1 -npfact 1 -mask ${FITSFNAME}_rfifind.mask -pdstep 2 -dmstep 1 -nopdsearch -fine -noscales -nooffsets -runavg -pd 0 ${FITSFNAME}.fits

#prepfold -nsub 288 -noxwin -ndmfact 1 -p ${PERIOD} -dm ${DM} -o ${OUTDIR}/${OBSID}/manual_fold_norunavg_${STAMP} -n 100 -pstep 1 -npfact 1 -mask ${FITSFNAME}_rfifind.mask -pdstep 2 -dmstep 1 -nopdsearch -fine -noscales -nooffsets -runavg -pd 0 ${FITSFNAME}.fits

#Creating the time series

if [ $(echo "${DM} < 40.52" | bc) -eq 1 ];then
    nmoutval=7392000
elif [[ $(echo "${DM} >= 40.52" | bc) -eq 1 && $(echo "${DM} < 116.42" | bc) -eq 1 ]];then
    nmoutval=3696000
else
    nmoutval=1848000
fi

prepsubband -numout ${nmoutval} -downsamp ${DS} -lodm ${DM} -dmstep 0.01 -numdms 1 -noclip -nsub 288 -mask ${RAWDIR}/${OBSID}/${beamtype}/SAP${SAP}/BEAM${BEAM}/${OBSID}_SAP${SAP}_BEAM${BEAM}_spdspsr_rfifind.mask -runavg -noscales -noweights -nooffsets -o ${OUTDIR}/${OBSID}/manual_fold_DM${DM}_downsamp_${DS}_${STAMP} ${RAWDIR}/${OBSID}/${beamtype}/SAP${SAP}/BEAM${BEAM}/${FITSFNAME}.fil

prepsubband -numout ${nmoutval} -downsamp ${DS} -lodm ${DM} -dmstep 0.01 -numdms 1 -noclip -nsub 288 -mask ${RAWDIR}/${OBSID}/${beamtype}/SAP${SAP}/BEAM${BEAM}/${OBSID}_SAP${SAP}_BEAM${BEAM}_spdspsr_rfifind.mask -noscales -noweights -nooffsets -o ${OUTDIR}/${OBSID}/manual_fold_DM${DM}_downsamp_${DS}_norunavg_${STAMP} ${RAWDIR}/${OBSID}/${beamtype}/SAP${SAP}/BEAM${BEAM}/${FITSFNAME}.fil
