#!/bin/bash


if [ $# -ne 1 ] ;then
    echo "This script reads the file with the positive candidates creted by the classifier, and copies the scrunched pfd files along with their png in the target location"
    echo ""
    echo "Usage: sh predictions.sh <obs id>"
    echo "i.e., sh predictions.sh L1212121"
    exit
fi

#obs ID
obsid="$1"

#input file location 
predfile="/projects/0/lotaas2/data/predictions/predictions_$obsid.txt"

#remove the file extensions from the list in the input file
predfilenames=`sed 's/.36scrunch//g' $predfile`

#output directory
outdir="/projects/0/lotaas2/data/pcands"

#create the observation directory
mkdir $outdir/$obsid

for file in $predfilenames;do
    cp $file.36scrunch  $outdir/$obsid
    convert $file.ps -rotate 90 $file.png
    mv $file.png $outdir/$obsid
done

cp $predfile $outdir/$obsid

#create the tarball for viewing
#mkdir ${outdir}/${obsid}_png

#cp $outdir/$obsid/*.png ${outdir}/${obsid}_png
#cd ${outdir}
#tar -cf ${obsid}_png.tar.gz ${obsid}_png
#rm -rf ${obsid}_png
#mv ${obsid}_png.tar.gz /projects/0/lotaas/data/out/check