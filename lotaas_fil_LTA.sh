#!/bin/bash -x
#SBATCH -p normal
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=sotiris.sanidas@gmail.com
#SBATCH --constraint=haswell

## if [ -z "$1" ] || [ "$1" = "-h" ]
##         then    
##         echo "Processes nbeams for an LOTAAS observation"
##         echo "USAGE: lotaas.sh TBD"
##         exit
## fi

PSRHOME=/home/joerivl


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

OUTDIR=/projects/0/lotaas2/data/out/$1/SAP$2
time {
   
# Set up some variables

    obsfile=$1    # e.g., L155450
    subarray=$2   # e.g., 1  for SAP1
    beam=$3       # e.g., 10 for BEAM10

    #set a dependent job to archive the scrunched pfd files
    #sbatch --dependency=afterany:$SLURM_JOB_ID /home/cooper/LOTAAS-Scripts/archive_pfd.sh $obsfile $subarray $beam
    
    node=`hostname`    
#RFI mask created directly from the .fil files - no master record
#    zapfile=frequencies4.zaplist
#    masterrfifile=/home/sanidas/LOTAAS-Scripts/master.rfi
    
    # Download file to local disk
    echo 'Copying data'
    
    # New location:
    mkdir BEAM${beam}; cd BEAM${beam}    
    
    stokes="stokes"
    
    time 
        
    if [[ $beam == 12 ]]
    then
        stokes="incoherentstokes"
    fi

    python $LOTAAS_PL/psrfits2fil.py --noscales --noweights --nooffsets /projects/0/lotaas2/data/raw/${obsfile}_red/$stokes/SAP${subarray}/BEAM${beam}/${obsfile}_SAP${subarray}_BEAM${beam}.fits -o ${obsfile}_SAP${subarray}_BEAM${beam}.fil

    
    # make new rfifind 
    rfifind -blocks 10 -noweights -noscales -nooffsets -o ${obsfile}_SAP${subarray}_BEAM${beam} *.fil

    cp *.mask $WORKDIR
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
	echo "outname is $outname"
	
	time srun mpiprepsubband \
	    -numout $nout -downsamp $ds -numdms $[$numberOFdms +1] -dmstep $dminc \
	    -noclip -nsub 288 -lodm $lodm -mask *.mask -runavg -noscales -noweights -nooffsets -o $outname_mpi *.fil


	date

        # single-pulse searches
        date
        for file in *.dat; do
                echo "./single_pulse_search_lotaas.py -m 0.105 -p $file " >> process_sp
        done
        date
        ./parallel.sh process_sp
        date
        time waitForFinish '[s]'ingle
        # create unique 'singlepulse' ascii file that will hold results for this iteration
        curdir=`pwd`
        spfile=`mktemp --tmpdir=$curdir singlepulse.XXXXXX`
        # now we combine .singlepulse files in one file
        echo "Concatenating .singlepulse to $spfile: cat *.singlepulse > $spfile"
        cat *.singlepulse > $spfile
        rm *.singlepulse

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
    echo "cat singlepulse.XXXXXX > ${obsfile}_SAP${subarray}_BEAM${beam}.singlepulse"
    cat singlepulse.* | head -n 1 > ${obsfile}_SAP${subarray}_BEAM${beam}.singlepulse
    cat singlepulse.* | grep -v \# >> ${obsfile}_SAP${subarray}_BEAM${beam}.singlepulse
    # making two separate master lists with widths <=20 and >=30
    cp ${obsfile}_SAP${subarray}_BEAM${beam}_rfifind.inf ${obsfile}_SAP${subarray}_BEAM${beam}.inf
    rm singlepulse.??????
    # making single-pulse plots
    echo "Making single-pulse plots for ${obsfile}_SAP${subarray}_BEAM${beam}..."
    dms=(0.0  9.0  19.0 20.0 100.0  300.0) # start DM for different ranges
    dme=(10.0  20.0 30.0 110.0  310.0  1000000.0) # end DM for different ranges
    stti=(0 600 1200 1800 2400 3000) # start times for different time intervals
    enti=(600 1200 1800 2400 3000 3600) # end times for different time intervals
    ws=(0.0  1.0  9.0 40.0 90.0) # start widths for different ranges
    we=(2.0  10.0 50.0 100.0 1000000.0) # end widths for different ranges
    for (( ii=0; ii<${#dms[@]}; ii++ )); do
        lowdm=`printf "%.2f" ${dms[$ii]}`
        highdm=`printf "%.2f" ${dme[$ii]}`
        echo "./single_pulse_search_lotaas.py -t 7.0 --dms ${dms[$ii]} --dme ${dme[$ii]} ${obsfile}_SAP${subarray}_BEAM${beam}.singlepulse" >> process_spplot
        # will also do several plots for different width ranges
        for (( kk=0; kk<${#we[@]}; kk++ )); do
         lowwth=`printf "%.2f" ${ws[$kk]}`
         highwth=`printf "%.2f" ${we[$kk]}`
         echo "./single_pulse_search_lotaas.py -t 7.0 --dms ${dms[$ii]} --dme ${dme[$ii]} --wmin ${ws[$kk]} --wmax ${we[$kk]} ${obsfile}_SAP${subarray}_BEAM${beam}.singlepulse" >> process_spplot
        done
        # will also do several plots for different time ranges
        for (( jj=0; jj<${#stti[@]}; jj++ )); do
         lowt=`printf "%.0f" ${stti[$jj]}`
         hight=`printf "%.0f" ${enti[$jj]}`
         echo "./single_pulse_search_lotaas.py -t 7.0 --dms ${dms[$ii]} --dme ${dme[$ii]} -s ${stti[$jj]} -e ${enti[$jj]} ${obsfile}_SAP${subarray}_BEAM${beam}.singlepulse" >> process_spplot
        done
    done
    date
    ./parallel.sh process_spplot
    date
    time waitForFinish '[s]'ingle
    for ps in *_singlepulse.ps; do
     base=`basename $ps .ps`
     echo "convert $ps $base.png" >> process_sppng
    done
    date
    ./parallel.sh process_sppng
    date
    time waitForFinish '[c]'onver
    tar cvfz ${obsfile}_SAP${subarray}_BEAM${beam}_singlepulse.tgz ${obsfile}_SAP${subarray}_BEAM${beam}.inf *.singlepulse *_singlepulse.ps *_singlepulse.png
    rm -f ${obsfile}_SAP${subarray}_BEAM${beam}.inf *.singlepulse *_singlepulse.ps *_singlepulse.png


    
    mkdir -p $WORKDIR/BEAM$beam"_sift"/sp
    # copying singlepulse tarball here
    cp ${obsfile}_SAP${subarray}_BEAM${beam}_singlepulse.tgz $WORKDIR/BEAM$beam"_sift"/sp

    cd $WORKDIR
    echo $PWD
   
            
    echo 'Python sifting and folding:'
    time python $PSRHOME/LOTAAS-Scripts/n_beam_sift_fold.py --ncores=24 --i_mixed=BEAM$beam --o=BEAM$beam"_sift" --i=BEAM$beam --filetype=FIL --masks=`ls BEAM$beam/*.mask`

    echo 'Memory usage at the end of a run:'
#    top -b|head -n 20
    df -h /dev/shm

    mkdir -p $outname_mpi
    mv BEAM$beam BEAM$beam"_sift" $outname_mpi 
    mv $WORKDIR/*.mask $outname_mpi
    cd $outname_mpi/BEAM$beam"_sift"/
   
    mkdir -p ${OUTDIR}

    #Scrunch pfd
    
    for cand in `find /dev/shm/tmp*/L*/*_sift/ -type f -name "*.pfd"`
    do
	echo "/home/joerivl/src/presto/bin/modify_pfd $cand 8" >> process_modify
    done
    /home/sanidas/LOTAAS-Scripts/parallel.sh process_modify
    date
    time waitForFinish '[m]'odify
    rm process_modify

    # Don't need to make a tar ball as this is done in archive script
    #tar -rf ${outname_mpi}_scrunched.tar `find FOLDS/ -type f -name "*.36scrunch"`
    #echo "${OUTDIR}/${outname_mpi}/BEAM${beam}_sift/${outname_mpi}_scrunched.tar" >> /home/cooper/LOTAAS-Scripts/scrunch2archive.lis
        
    find /dev/shm/tmp*/L*/*_sift/ -type f -size +1000k -name "*.pfd" -delete
    

    cd ../BEAM$beam/
    rm *.fil
    tar -cf ${outname_mpi}.tar ${outname_mpi}_DM* /home/sanidas/LOTAAS-Scripts/lotaas_fil_LTA.sh 
    
    #mv $outname_mpi.tar ../
    rm -f `ls * | grep -v tar`
    rm process*

    echo 'Total time:'
}

# copy data back
set -x 
mkdir -p ${OUTDIR}
echo "What is in wor
     kdir ================== "
#find ${WORKDIR}

cp -r ${WORKDIR}/* ${OUTDIR}
cd $OUTDIR

#Score generator Classifier 
python $LOTAAS_PL/scores_robLyon/PulsarProcessingScripts-master/src/CandidateScoreGenerators/ScoreGenerator.py -c $OUTDIR/$outname_mpi/ --pfd --dmprof --arff -o $OUTDIR/$outname_mpi/dmprof_$beam.arff 


rm -Rf ${WORKDIR}
  
