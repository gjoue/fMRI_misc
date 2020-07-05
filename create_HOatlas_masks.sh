FSLDIR=/common/apps/fsl

## stuff to get fsl working cos not setup on the server
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH

thr=25;
res=2;

#dirOut=~/distorcorr/atlasMasks
dirOut=/projects/crunchie/estropharm3/masks_HOatlas
dirAtlases=$FSLDIR/data/atlases
dirHO=${dirAtlases}/HarvardOxford
cort=${dirHO}/HarvardOxford-cort-maxprob-thr${thr}-${res}mm.nii.gz
subcort=${dirHO}/HarvardOxford-sub-maxprob-thr${thr}-${res}mm.nii.gz

file_cortLab=/projects/crunchie/TOOLBOXES/ATLASES/HOatlasMasks/HarvardOxford-Cortical_ind2abbr.txt
file_subcortLab=/projects/crunchie/TOOLBOXES/ATLASES/HOatlasMasks/HarvardOxford-Subcortical_ind2abbr.txt

## copied the Harvard-Oxford cortical atlas ROI IDs from
## https://scalablebrainatlas.incf.org/services/labelmapper.php?template=HOA06

readarray lab_cort < $file_cortLab
readarray lab_subcort < $file_subcortLab

for l in ${!lab_cort[*]}
do
    roientry=${lab_cort[$l]}
    echo .........Creating mask from $roientry
    ind=$(echo $roientry | awk '{print $1}')
    roi=$(echo $roientry | awk '{print $2}')
    echo $roi has index $ind

    fslmaths $cort -thr $ind -uthr $ind $dirOut/HOcort_${roi}_thr${thr} 
done

for l in ${!lab_subcort[*]}
do
    roientry=${lab_subcort[$l]}
    echo .........Creating mask from $roientry
    ind=$(echo $roientry | awk '{print $1}')
    roi=$(echo $roientry | awk '{print $2}')
    echo $roi has index $ind

    fslmaths $subcort -thr $ind -uthr $ind $dirOut/HOsubcort_${roi}_thr${thr} 
done



##................... ATTEMPT TO READ THE XML THAT COMES WITH HARVD-OXFD ATLAS.........
## need to add 1 to the index values in 
## /common/apps/fsl/data/atlases/HarvardOxford-Subcortical.xml 
## 


#subcortLab="${dirHO}/HarvardOxford-Subcortical.xml"
#cortLab="${dirAtlases}/HarvardOxford-Cortical.xml"

## bash will separate items in array with spaces -- unset this because there are spaces in the labels
#IFS=$'\t\n' # set (Internal Field Separator) variable so that it splits fields by tab and newline not spaces

#labs=( $(xmlstarlet sel -T -t -m //data/label  -c . -o ";" -v @index -n "$cortLab"))
#echo Found  ${#labs[@]} label

## https://www.joyofdata.de/blog/transforming-xml-document-into-csv-using-xmlstarlet/
## count number of labels in file and print new line -n
## -t : the following parameters are part of the output template
## -v : output the value of an xpath expression
## -m : iterate over all nodes that match the provided xpath expression
#  xmlstarlet sel -t -v "count (//label)" -n HarvardOxford-Subcortical.xml

#for i in ${!labs[*]}
#do
#	  echo "$i" "${labs[$i]}"
#	    # instead of echo use the values to send emails, etc
#    done

