#!/bin/bash -x
#SBATCH -p normal
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=sally.cooper@postgrad.manchester.ac.uk
#SBATCH --constraint=haswell

## if [ -z "$1" ] || [ "$1" = "-h" ]
##         then    
##         echo "Processes nbeams for an LOTAAS observation"
##         echo "USAGE: lotaas.sh TBD"
##         exit
## fi

PSRHOME=/home/joerivl
echo "START TIME"
echo `date`
################################################################################
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
	while [ $job -gt 0 ] ; do sleep 1; top -b | head -n 40; job=`ps ax | grep ${STRING} | wc -l`; done	
    done
}
################################################################################


#WORKDIR=$TMPDIR
rm -rf /dev/shm/tmp/*
WORKDIR=`mktemp -d --tmpdir=/dev/shm`
cd $WORKDIR
pwd
echo "START=" `date`

OUTDIR=/projects/0/lotaas2/data/out/confirmations/$1/SAP$2
time {
   
# Set up some variables

    obsfile=$1    # e.g., L155450
    subarray=$2   # e.g., 1  for SAP1
    beam=$3       # e.g., 10 for BEAM10
#   beam=$(( ( RANDOM % 74 ) )) # For testing purposes on large N Beams
    
#    node=`hostname`    
#    zapfile=frequencies4.zaplist
#    masterrfifile=/home/sanidas/LOTAAS-Scripts/master.rfi
    
    # Download file to local disk
    echo 'Copying data'
    
    # New location:
    mkdir BEAM${beam}; cd BEAM${beam}    
    
    stokes="stokes"
    
    time 
        
    if [[ $beam == 0 ]]
    then
	stokes="incoherentstokes"
    fi

    python $LOTAAS_PL/psrfits2fil.py --noscales --noweights --nooffsets /projects/0/lotaas2/data/raw/confirmations/${obsfile}_red/$stokes/SAP${subarray}/BEAM${beam}/${obsfile}_SAP${subarray}_BEAM${beam}.fits -o ${obsfile}_SAP${subarray}_BEAM${beam}.fil

    
    # make new rfifind 
    rfifind -blocks 10 -noweights -noscales -nooffsets -o ${obsfile}_SAP${subarray}_BEAM${beam} *.fil

    ls 

        
    # processing using shell loop ...
    dminc=0.01
    dm1inc=0.05
    dm2inc=0.1
    dmulim1s=3976
    dmulim2s=6000
    #numberOFdms=252
    num=0

    cp $LOTAAS_PL/parallel.sh . 
    cp $LOTAAS_PL/realfft.sh . 
    cp $LOTAAS_PL/rednoise.sh .
    cp $LOTAAS_PL/single_pulse_search_lotaas.py .
    
    #for num in $(seq 0 $numberOFdms 9600); do 
    while [ $num -le  9600 ]
    do
	# The various steps
	# Step 1
	if [ $num -lt $dmulim1s ] 
	    then
	    numberOFdms=252
	    lodm=`echo "$num*$dminc" | bc`
	    dmulim1=`echo "$num + $numberOFdms" | bc`
	    dm1=`echo "$dmulim1*$dminc" | bc`
 	    ds=1
 	    nout=7392000 # output from choose_N
	    let num=num+253
	# Step 2
	elif [ $num -gt $dmulim1s ] && [ $num -lt $dmulim2s ]
	    then
	    numberOFdms=505
	    lodm=`echo "$dm1 + ($num-$dmulim1)*$dm1inc" | bc`
	    dmulim2=`echo "$num + $numberOFdms" | bc`
	    dm2=`echo "$dm1 + ($dmulim2-$dmulim1)*$dm1inc" | bc`
	    dminc=$dm1inc
	    ds=2
	    nout=3696000
	    let num=num+506
	# Step 3
	elif [ $num -gt $dmulim2s ]
	    then
	    numberOFdms=1011
	    lodm=`echo "$dm2 + ($num-$dmulim2)*$dm2inc" | bc`
	    dminc=$dm2inc
	    ds=4
	    nout=1848000
	    let num=num+1012
	else
	    # we shouldn't get here
	    echo "Error -- we shouldn't get here"
	fi	

	rm -f process_sp process_spplot process_sppng process_fft process_red process_accel*
   
	outname_mpi=`ls *.fil | awk -F\. '{print $1}'`
	
	time srun mpiprepsubband \
	    -numout $nout -downsamp $ds -numdms $[$numberOFdms +1] -dmstep $dminc \
	    -noclip -nsub 288 -lodm $lodm -mask *.mask -runavg -noscales -noweights -nooffsets -o $outname_mpi *.fil


	date

        # single-pulse searches

	count=0   # the number of commands per spawn (i.e., 92 DM/24 cores => 4)
        bcount=0
	date
	for file in *.dat; do
	    if [ $count -lt 3 ] && [ $bcount -lt $numberOFdms ] # 4 commands per call, so count < 3
		then
#		echo -n "realfft -outdir . $file > $file.log; " >> process_fft
		echo -n "./realfft.sh $file; " >> process_fft
		let count+=1
	    else
#		echo    "realfft -outdir . $file > $file.log  " >> process_fft
		echo "./realfft.sh $file " >> process_fft
		count=0
	    fi
            let bcount+=1
	done
	date

	./parallel.sh process_fft
	date
	waitForFinish '[r]'ealff
	
	count=0
        bcount=0
	for file in *.fft; do
	    if [ $count -lt 3 ] && [ $bcount -lt $numberOFdms ]
		then
#		echo -n "rednoise $file >> $file.log ;" >> process_red
		echo -n "./rednoise.sh $file ;" >> process_red
		let count+=1
	    else
#		echo "rednoise $file >> $file.log" >> process_red
		echo "./rednoise.sh $file" >> process_red
		count=0
	    fi
	    let bcount+=1
	done
	
	./parallel.sh process_red; 
	date
	waitForFinish '[r]'edno
	
        #rm -f `ls *.fft | grep -v red`	
	for file in `ls *red.fft | awk -F"_red.fft" '{print $1}'`; do
	    file_long=$file"_red.fft"
	    echo "mv  $file_long $file.fft" >> process_mv
	done
	source process_mv
	rm process_mv
	
	count=0
	bcount=0
        ccount=0
	for file in *.fft; do
	    if [ $count -lt 3 ] && [ $bcount -lt $numberOFdms ]
                then
		echo -n "accelsearch -flo 1. -sigma 2 -fhi 1000. -zmax 0 -numharm 16 -harmpolish $file > $file.log;" >> process_accel_$ccount
		let count+=1
	    else
		echo    "accelsearch -flo 1. -sigma 2 -fhi 1000. -zmax 0 -numharm 16 -harmpolish $file > $file.log " >> process_accel_$ccount
		count=0
		echo "source ./process_accel_$ccount" >> process_accel
		let ccount+=1
	    fi
	    let bcount+=1
	done

	./parallel.sh process_accel
	echo 'time waitForFinish accelsearch:'
	time waitForFinish '[a]'ccel 


	date
	rm -f $outname_mpi*.dat $outname_mpi*.fft
##	mkdir -p save; mv $outname_mpi*.dat $outname_mpi*.fft save/

    done #The big loop over the DM steps

    # merging all singlepulse.XXXXXX files together

    cd $WORKDIR
    echo $PWD
   
    
    
    
        
    echo 'Python sifting and folding:'
    time python $PSRHOME/LOTAAS-Scripts/n_beam_sift_fold.py --ncores=24 --i_mixed=BEAM$beam --o=BEAM$beam"_sift" --i=BEAM$beam --filetype=FIL --masks=`ls BEAM$beam/*.mask`

    echo 'Memory usage at the end of a run:'
#    top -b|head -n 20
    df -h /dev/shm

    

    mkdir -p $outname_mpi
    mv BEAM$beam BEAM$beam"_sift" $outname_mpi 
    
    cd $outname_mpi/BEAM$beam"_sift"/
   
    mkdir -p ${OUTDIR}

    #Score generator Classifier
    python $LOTAAS_PL/scores_robLyon/PulsarProcessingScripts-master/src/CandidateScoreGenerators/ScoreGenerator.py -c . --pfd --dmprof --arff -o $OUTDIR/dmprof_$beam.arff 

    cp  $OUTDIR/dmprof_$beam.arff $OUTDIR/L*_SAP0_BEAM$beam/
    #Scrunch pfd
    
    for cand in `find /dev/shm/tmp*/L*/*_sift/ -type f -size +1000k -name "*.pfd"`
    do
	echo "/home/joerivl/src/presto/bin/modify_pfd $cand 8" >> process_modify
    done
    $LOTAAS_PL/parallel.sh process_modify
    date
    time waitForFinish '[m]'odify
    rm process_modify

#    tar -rf ${outname_mpi}_${beam}scrunched.tar `find /dev/shm/tmp*/L*/*_sift/ -type f -name "*.36scrunch"`
    
#    mkdir /archive/pulsar/lotaas/scrunchedCands/${outname_mpi}
#    cp ${outname_mpi}_${beam}scrunched.tar /archive/pulsar/lotaas/scrunchedCands/${obsid}/

    find /dev/shm/tmp*/L*/*_sift/ -type f -size +1000k -name "*.pfd" -delete
#    rm ${outname_mpi}_${beam}scrunched.tar

    cd ../BEAM$beam/
    rm *.fil
    tar -cf $outname_mpi.tar $outname_mpi_DM*
    rm -f `ls * | grep -v tar`
    rm process*

    echo 'Total time:'
}

# copy data back
set -x 
mkdir -p ${OUTDIR}
echo "What is in wor
     kdir ================== "
find ${WORKDIR}

cp -r ${WORKDIR}/* ${OUTDIR}
cd

rm -Rf ${WORKDIR}
  
echo "END TIME"
echo `date`