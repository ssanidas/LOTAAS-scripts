#!/bin/bash

# Call the scripts to produce diagnostic plots of specified candidate

# USAGE: sh candidates_heatmap.sh candID, fileID, where candID is the first column of LSP_candidates spreadsheet and fileID the first column of LTA spreadsheet
# EXAMPLE: sh candidates_heatmap.sh L101010_1_1_10 L101012

if [ $# -ne 2 ] ;then
    echo USAGE: sh candidates_heatmap.sh candID, fileID, where candID is the first column of LSP_candidates spreadsheet and fileID the first column of LTA spreadsheet
    echo EXAMPLE: sh candidates_heatmap.sh L101010_1_1_10 L101012
    exit
fi

candID=$1
fileID=$2

arr=(${candID//_/ })
OBS=${arr[0]}
SAP=${arr[1]}
BEAM=${arr[2]}
CAND=${arr[3]}


sh candidates_heatmap.sh $1
sh dedispersed_plot.sh $1 $2



