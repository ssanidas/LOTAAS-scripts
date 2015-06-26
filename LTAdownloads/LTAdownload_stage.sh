#! /bin/bash
#SBATCH -p staging
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-user=sotiris.sanidas@gmail.com

# Verbose mode
set -x

####################################
# wait to parallel script to finish
function waitForFinish()
{
    local STRING;

    STRING=$1;

    # wait until jobs have started
    sleep 1

    # check, twice, whether all done
    for i in 1 2 ; do
        job=99
        while [ $job -gt 0 ] ; do sleep 1; top -b | head -n 40; job=`ps ax | grep ${STRING} | wc -l`; 
done      
    done
}
####################################

#probably obsolete? The proxy file is also in the same directory with the script btw
#export_X509_USER_PROXY=/home/hessels/proxydir/proxyfile



# Read command line
IDS="$@"

# Set up Grid credentials
module load globus
echo P4ls@rgrid | voms-proxy-init -pwstdin -valid 168:00 -dont-verify-ac -voms lofar:/lofar/user
GSIDIR=gsiftp://gridftp.grid.sara.nl
BASEDIR=/projects/0/lotaas/data/raw/

for ID in $IDS ; do 
    i=0
    # Make dir
    LOCALDIR=${BASEDIR}/${ID}_red/
    mkdir -p ${LOCALDIR}
    cd ${LOCALDIR}

    rm process_globus
    skip=0
    while read line
    do

	#if (($skip >= 194 ))
	#then
	FILE=`echo $line | awk -F":8443/" '{print $2}' | awk '{print $1}'` 
        OUTFILE=`echo $line | awk -F/ '{print $15}'`
	
	echo "globus-url-copy -rst ${GSIDIR}${FILE} - | tar Bxp" >> process_globus
	#fi
	#let skip+=1
	
    done  < /home/sanidas/LTAdownloads/${ID}srm.txt


    /home/sanidas/LTAdownloads/parallelglobus.sh process_globus
    
    date
    waitForFinish '[g]'lobus


done
