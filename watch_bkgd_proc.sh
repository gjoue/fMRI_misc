## INPUT: associative array with PIDs as index and subjID as value
## WHAT IT DOES: prints messages
## for some reason can't get it to print all elements in the array??
##
## to kill all bkgd processes even when ctl+c on script (won't otherwise):
## trap "kill 0" EXIT
##
## EXAMPLE USAGE:
## for subj in "${doSubjs[@]}"
## do
##    ## do sth
##    pid=$!
##    pids[$pid]=$subj
##    echo "queued $subj - $pid"
##    ## waiting every iteration doesn't max use of multicores so limit job #s
##    joblist=($(jobs -p))
##    while (( ${#joblist[*]} >= 10 ))
##    do
##     	  echo ....maxing 10 proc, so sleeping 30s before pushing more onto queue
##     	  echo
##     	  sleep 30s
##     	  joblist=($(jobs -p))
##     done
##     $* & # $* expands to a single arg w/all elem sep'd b spaces
##     date
## done
## ## wait for all bkgd jobs to finish
## watch_bkgd_proc "pids"


watch_bkgd_proc () {
    #    pids=("$@")
    pids=$(declare -p "$1") # show associative array definition

    NPIDS=0
    NPIDS_1=0

    ## process ID = 0 in /dev/null but do not try to unset
    
    
    while [ -n "${pids[*]}" ]; do
	sleep 30s
	
	for pid in "${!pids[@]}"; do
	    if ! ps "$pid" >/dev/null ; then
#		if  $pid -ne 0 ; then
		    unset pids[$pid]
		    NPIDS="${#pids[@]}"
		    echo "unset: $pid"
#		fi
	    fi
	done
	if [ -z $NPIDS ]; then
	    break
	fi
	echo

	echo "Still waiting for $NPIDS processes to finish ($(date +%T)) ...."
	
	if [ $NPIDS -lt 5 ]  && [ $NPIDS -ne $NPIDS_1 ]; then
	    paste <(printf "PID %s\n" "${!pids[@]}") <(printf "(%s)\n" "${pids[@]}")
	    NPIDS_1=$NPIDS
	fi

    done

    }
