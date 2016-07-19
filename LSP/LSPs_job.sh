#!/bin/bash -x
#SBATCH -p normal 
#SBATCH -t 10:00:00
#SBATCH -n 24
#SBATCH --mail-user=danielemichilli@gmail.com


##########################################################
#
# LOTAAS Single Pulse search job launcher
#
# written by Daniele Michilli
#
# Not written for direct use, use LSPs_launcher.sh instead
#
##########################################################

python $LOTAAS_PL/../LSP/LSPs/bin/SPclean $1 -pl #-conf -folder $2
#python /home/danielem/scripts/alerts.py $2$1/sp

exit 0
