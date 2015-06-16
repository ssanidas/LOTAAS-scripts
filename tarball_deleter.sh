#!/bin/bash                                                                                                 
#This script picks up all the observations that have an sp tarball and deletes them
#USE WITH CAUTION!!!!!!!!!!!!!!!!!!!!
#ONLY WHEN EVERYTHING IS ON TAPE+VERIFIED

datadir="/projects/0/lotaas/data/out/tarballs"

cd $datadir

delobs=`ls |grep sp|awk -F_ '{print $1}'`


for i in $delobs; do
    echo "$i" >>/projects/0/lotaas/data/out/removed_tarballs
    echo "deleting $i""_*.tar"
    rm $i"_"*.tar
done