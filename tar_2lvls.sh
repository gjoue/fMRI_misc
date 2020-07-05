#!/bin/bash
## dir_lev1_a
##     dir_lev2_a
##     dir_lev2_b
## dir_lev1_b
##     dir_lev2_a
##     dir_lev2_c
##
## tar one level of subfolders in each folder as a separate tar ball
## tar files are labeled by the top level and bottom level folder names

##...... vars to modify ...............
baseDir='/projects/crunchie/bnst' # up to path that you don't want in the tar ball when it unpacks
sourceTar='DATA_MRI/DATA_bnst2018'
targetTar='data_orig/data_orig_bnst2018'

cd $baseDir
pwd


[[ -d $targetTar ]] || mkdir $targetTar

#declare -a subList=('sub047' 'sub048' 'sub049' 'sub050' 'sub051');

## in this example, each folder has the format sub%03d, so we can specify a range of numbers of the
## matching folders we want to go through
start=34;#31;#27;
end=40;#33;#30;

# for i in {0..2} {4..6} {7..11..2}; do echo $i; done
for i in $(seq $start $end)  
#for ss in "${subList[@]}"  # if you have specified a list of subfolders to go thru specifically
#for ss in `ls ${sourceTar}/indiv | grep sub`  # if you want to go through all subfolders
do
    ss=$(printf "sub%02d" "$i")
    subSrc="${sourceTar}/${ss}"

    echo ............. tarring subj $ss ............
    echo subsrc=$subSrc
    
    for ww in `ls ${subSrc}`
    do
	subfolder=${subSrc}/${ww} 

	ff="${ss}.${ww}.tar.gz"

	echo Checking archive $ff -- if does not exist then creating 
	    
	[[ -f ${targetTar}/${ff} ]] || tar cvfz $targetTar/$ff ${subfolder} > ${targetTar}/${ss}.${ww}.log
    done
done
