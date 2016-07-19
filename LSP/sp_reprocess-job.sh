#!/bin/bash -x

work_dir="/projects/0/lotaas2/data/out/sp"

cd $work_dir

for idL in "$@"; do
  ./extractor_for_sp.scr obs $idL

  sh ~/scripts/LSPs_launcher.sh $idL
  
  ./archiver.scr obs $idL sp


done








