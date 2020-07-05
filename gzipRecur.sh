## call with 1 input arg
## e.g. /projects/estropharm3/data/MRI/indiv
# s0=1
# sf=53


# while [ $s0 -le $sf ]
# do
#     subj="$(printf '%03d' $s0)"
#     echo sub=$subj
#     swd="${1}/sub${subj}"

#     for t in $swd/func/*
#     do
# 	for r in ${t}/run*
# 	do
# 	    echo zipping up files in $r;
# 	    find $r -name fPRISMA*.nii -exec gzip {} \;
# 	    find $r -name afPRISMA*.nii -exec gzip {} \;
# 	done
#     done
#     s0=$((s0+1))
# done
		      




for f in /projects/crunchie/dinu/catlearn1/fMRI_data/p06/*
do
    what=${f#/*p06/}
    echo $f, $what
    find $f -type d -name "original_fPRISMAs" -exec tar --exclude='before' --exclude='after' --exclude='dummy' -rvf ${what}.tar {} + > catlearn1_origPRISMA_${what}.out
done


#find /projects/crunchie/dinu/catlearn1/fMRI_data/p06/ -type d -name "run*" -exec tar -cvfz test.tar.gz {} \+
    
for f in /projects/crunchie/dinu/catlearn1/fMRI_data/p06/*
do
    what=${f#/*p06/}
    echo $f, $what
    find $f -type d -name "original_fPRISMAs" -exec tar --exclude='before' --exclude='after' --exclude='dummy' -rvf ${what}.tar {} + > catlearn1_origPRISMA_${what}.out
done
