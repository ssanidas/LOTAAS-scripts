#!/bin/bash

#SP Candidate Periodicity scanner
#A simple script that gets a particular beam from the archive and 
#prepares all the presto periodicity candidates png files

if [ $# -ne 4 ] ;then
    echo "The script needs as arguements the obsid, sap, beam and DM of the candidate to investigate"
    echo "i.e., sh spcp_scan.sh L121212 1 45 5"
    exit
fi


OBS=$1 #obs id
SAP=$2 #sap
BEAM=$3 #beam
DM=$4 #candidate's DM

OUTDIR=/projects/0/lotaas2/data/out/check
WRKDIR=/projects/0/lotaas2/data/out

#fetch the tarball
cd /archive/pulsar/lotaas/out_tarballs
echo "Tarballs:"
echo `dmls -l ${OBS}_SAP${SAP}*|awk '{print $9}'`
arcfile=`dmls -l ${OBS}_SAP${SAP}*|awk '{print $9}'`
dmget -a $arcfile

#extract the desired beam
mkdir ${WRKDIR}/${OBS}
tar -xf $arcfile -C ${WRKDIR}/${OBS} SAP${SAP}/${OBS}_SAP${SAP}_BEAM${BEAM}

#create the pngs

cd ${WRKDIR}/${OBS}/SAP${SAP}/${OBS}_SAP${SAP}_BEAM${BEAM}/BEAM${BEAM}_sift/FOLDS

mkdir ${OUTDIR}/${OBS}
 
find . -name "*.ps" >filelist_tmp
dmmin=$((DM-2))
dmmax=$((DM+2))
if [ -f filelist_clean_tmp ];then 
    rm filelist_clean_tmp
    touch filelist_clean_tmp
else
    touch filelist_clean_tmp
fi

for i in `seq ${dmmin} ${dmmax}`;do    
    cat filelist_tmp|grep -e "DM${i}\." >>filelist_clean_tmp
done


for i in `cat filelist_clean_tmp`;do
   fname=`echo "$i"|awk -F/ '{print $3}'  |sed 's/.ps//g'`
   convert $i -rotate 90 ${OUTDIR}/${OBS}/${fname}.png
done

rm filelist_tmp 
rm filelist_clean_tmp

cd ${OUTDIR}
tar -cf ${OBS}.tar.gz ${OBS}
rm -rf ${OBS}

rm -rf ${WRKDIR}/${OBS}