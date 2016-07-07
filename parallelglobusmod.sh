#!/bin/bash

# Maximum number of jobs


while getopts ":c:" opt;do
    case $opt in
        c)
            MAXJOBS=$OPTARG
            if [ "$MAXJOBS" \> "8" ] || [ "$MAXJOBS" \< "1" ];then
                echo "Option $OPTARG not allowed. Must be an integer between 1 and 8."
                exit
            fi
            ;;
        \?)
            echo "Invalid Option: -$OPTARG" >&2
            echo "Currently allowed options: -t,c"
            exit
            ;;
        :)
            echo "Option -$OPTARG must have an argument"
            exit
            ;;
    esac
done
shift $((OPTIND-1))





# Get job ID from scrint
function jobidfromstring()
{
    local STRING;
    local RET;

    STRING=$1;
    RET="$(echo $STRING | sed 's/^[^0-9]*//' | sed 's/[^0-9].*$//')"

    echo $RET;
}

# Ready to spawn job
function clearToSpawn
{
    local JOBCOUNT="$(jobs -r | grep -c .)"
    if [ $JOBCOUNT -lt $MAXJOBS ] ; then
        echo 1;
        return 1;
    fi

    echo 0;
    return 0;
}

# Loop over commands
JOBLIST=""
cat $1 | while read line;
  do
  while [ `clearToSpawn` -ne 1 ] ; do
      sleep 0.5
  done

  # Run command
  echo "Running" $line
  echo $line | sh &

  LASTJOB=`jobidfromstring $(jobs %%)`
  JOBLIST="$JOBLIST $LASTJOB"
done
# Wait for jobs to finish
for JOB in $JOBLIST ; do
    wait %$JOB
done
