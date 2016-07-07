#!/bin/bash

#A quick and dirty script to fix the coordinate issues in the early lotaas observations - this has to be executed in the same directory as the pointing

#Usage: sh coor_fixer.sh L121212

if [ $# -ne 1 ];then
    echo "A shamelessly quick n dirty script to fix the coordinate issues  in already downloaded LOTAAS pointings. "
    echo ""
    echo "Usage: sh coor_fixer.sh L121212"
    echo ""
    echo "Warning:It must be executed at the same directory as the pointing"
    exit
fi

pointing=$1
repository=/projects/0/lotaas/software/LOTAAS-Scripts
rawdir=/projects/0/lotaas/data/raw

cd "$pointing""_red"
echo "I am in "`pwd`

for i in {0..2};do
    cd stokes/SAP$i
    for j in {0..67};do
	cd BEAM$j
	python ${repository}/fix_fits_coordinates.py -v -d ${rawdir}/$pointing"_red" $pointing"_SAP"$i"_BEAM"$j.fits
	if [[ $BEAM == 6 ]];then
	    python ${repository}/fix_coords_IS_beams.py --meta ../../../incoherentstokes/SAP$i/BEAM$j/*.h5 --verbose *.fits
	fi
	cd ..
    done
    cd ../..
done