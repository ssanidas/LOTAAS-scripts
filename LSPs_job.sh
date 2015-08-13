#!/bin/bash -x
#SBATCH -p normal
#SBATCH -t 5:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1



##########################################################
#
# LOTAAS Single Pulse search job launcher
#
# written by Daniele Michilli
#
# Not written for direct use, use LSPs_launcher.sh instead
#
##########################################################

python LSPs/bin/SPclean $1 $2 -pl

exit 0
