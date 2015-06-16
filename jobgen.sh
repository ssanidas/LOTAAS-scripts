#!/bin/bash

#the jobs scripts directory
jobdir="/home/sanidas/jobdir"

if [ ! -d $jobdir ];then
    mkdir $jobdir
    else
    rm -r $jobdir
    mkdir $jobdir
fi

dirs="L182997" 

for dir in $dirs;do
    if [[ $dir == L* ]];then
	echo "#!/bin/bash" >>$wdir/cpconfjobs/job_$dir.scr
	echo "#SBATCH -p staging" >>$wdir/cpconfjobs/job_$dir.scr
	echo "#SBATCH -t 72:00:00" >>$wdir/cpconfjobs/job_$dir.scr
	echo "sh copy_confs.sh $dir" >>$wdir/cpconfjobs/job_$dir.scr
    fi
done


