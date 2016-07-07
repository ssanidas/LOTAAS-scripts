#!/bin/bash

#==================================================#
#Heatmap plot generator for processed confirmations#
#==================================================#

#Note 1
#This script has to be executed inside the directory of the observation

if [ $# -ne 1 ] && [ $# -ne 2 ] && [ $# -ne 3 ];then
    echo "This script produces the heatmap for a specific candidate on a processed observation confirmation. The script requires three user provided arguements. The first is the DM of the candidate. It will pick up the beams that it was seen, creating a heatmap based on the chi-squared values that were produced by the processing pipeline. The second arguement is the name of the candidate, as this will appear in the plots. The third arguement, is the obsid. The processed and raw data have to be in fixed locations (these are the standard cartesius locations;the user has to change these to make the script working at a different environment). The script will pick up any source with a DM +/-2 of the one provided by the user."
    echo ""
    echo "        i.e., sh heatmap_conf.sh 25 J1010+10 L121212"
    echo ""
    echo "The script has to be executed inside the directory of the processed follow-up observation. Good hunting!"
    exit
fi



#The obsID
obsid=$3
#obsid=`pwd|awk -F/ '{print $NF}'`
#fullobsid=`pwd`
fullobsid=/projects/0/lotaas2/data/out/confirmations/$obsid

echo $obsid
echo $fullobsid
#raw data directory
rawdata=/projects/0/lotaas2/data/raw/confirmations

#the psr name
psrname=$2

#the target DM
dm=$1

#min DM value
dmmin=$(($1-2))
echo "DM min:$dmmin"
#max DM value
dmmax=$(($1+2))
echo "DM max:$dmmax"

#remove chi-squared_scripted file if present
if [ -e $fullobsid/chi-squared_scripted.txt ];then
    rm $fullobsid/chi-squared_scripted.txt
fi

cd $fullobsid/SAP0

for i in {1..127};do #change to 1 and 127
    if [ -d $fullobsid/SAP0/$obsid"_SAP0_BEAM"$i ];then
	cd $fullobsid/SAP0/$obsid"_SAP0_BEAM"$i
	cd "BEAM"$i"_sift/FOLDS"
	for j in {0..23};do #change to 23 #searching inside the cpu-core directories
	    cd "CORE_"$j
	    filenames=`ls|grep bestprof`
	    for k in $filenames;do
		dmval=` echo "$k" |awk -F_ '{print substr($4,3)}'|awk -F. '{print $1}'`
		if [ "$dmval" -lt "$dmmax" ] && [ "$dmval" -gt "$dmmin" ];then
		    snr=`cat $k |grep "Reduced"|awk '{print $5}'`
		    echo "$k $snr $j $i" >>../templist_delete
		    else
		    echo "$k 0.1 $j $i" >>../templist_delete
                fi
	    done
	    cd ..
	done
	selection=`cat templist_delete|sort -n -k 2|tail -n1`
	filename=`echo "$selection"|awk '{print $1}'|sed 's/bestprof/ps/g'`
	corename=`echo "$selection"|awk '{print $3}'`
	beamname=`echo "$selection"|awk '{print $4}'`
	snr=`echo "$selection"|awk '{print $2}'`
	partname=`echo $filename|awk -F_ '{print "_"$4"_"$5"_"$6"_"$7"_"$8}'`
	imgname="file=$fullobsid/SAP0/$obsid""_SAP0_BEAM$beamname/BEAM$beamname""_sift/FOLDS/CORE_$corename/$psrname""_$obsid""_SAP0_BEAM$beamname""_PSR_$psrname$partname"


	echo "$imgname obs=CS_SAP0_BEAM$beamname""_$psrname S/N=$snr" >>$fullobsid/chi-squared_scripted.txt
	rm templist_delete
    else
	echo "Directory "$obsid"_SAP0_BEAM"$i" does not exist"
    fi
done


#create the heatmap plot
python $LOTAAS_PL/plot_LOFAR_TA_multibeam3.py -c $fullobsid/chi-squared_scripted.txt --target $psrname --sap 0 --parset $rawdata/$obsid"_red"/*.pseudo.parset --out_logscale $fullobsid/$obsid"_heatmap_log.png" --out_linscale $fullobsid/$obsid"_heatmap_lin.png"

convert -append $fullobsid/$obsid"_heatmap_log.png" $fullobsid/$obsid"_heatmap_lin.png" $fullobsid/TAheatmap_${obsid}_search.png

#copy the parset file
cp $rawdata/$obsid"_red"/*.pseudo.parset $fullobsid
cp $rawdata/$obsid"_red"/*.par $fullobsid