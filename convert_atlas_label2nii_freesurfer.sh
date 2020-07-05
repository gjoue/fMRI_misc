

## This shows the commands of how a binary mask was created from
## FreeSurfer labels

## A word on FreeSurfer files:
##
## fsaverage = avg of 40 Ss using spherical averaging (Fischl et al. 1999)
## fsaverage, fsaverage6, fsaverage5 = Resting State Cortical Parcellation
##  https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011
## "fsaverage" contains the high resolution version of the
## parcellation, while "fsaverage6" and "fsaverage5" contain lower
## resolution versions of the parcellation. The parcellations were
## computed in "fsaverage5" space and upsampled to "fsaverage6" and
## "fsaverage".
##
## for attempts to map between volumetric (e.g. MNI152) to surface
## (e.g. fsaverage) spaces:
## https://ww5.aievolution.com/hbm1701/index.cfm?do=abs.viewAbs&abs=1236
##
## fsaverage/mri.2mm
##   - contains volumes sampled into 2mm space for volume-based fMRI analysislr
## $FREESURFER_HOME/average/mni152.register.dat
##   - reg between fsaverage/mni305 subj (full 256^3, 1mm3) and 2mm mni152 sp
export FREESURFER_HOME="/common/apps/freesurfer-6.0.0_x86_64"
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR="."
FSLDIR="/common/apps/fsl"

#> FSFAST_HOME       /Applications/freesurfer/fsfast
#> FSF_OUTPUT_FORMAT nii.gz
#SUBJECTS_DIR="${FREESURFER_HOME}/subjects"
#> MNI_DIR           /Applications/freesurfer/mni



# ##..... transform labels into some space
# # -reg creates a registration matrix file
# tkregister2 –mov /Applications/freesurfer/subjects/${blindnum}/mri/rawavg.mgz \ #input
# 	    –noedit \           # suppresses GUI
# 	    –s ${blindnum} \    #input subj
# 	    –regheader \        # read hdr from input
# 	    –reg ./register.dat # creates a reg matrix

# ##....... convert fsaverage to MNI space

# mri_vol2vol --targ fsaverage.DKT40+aseg.nii \
# 	    --mov $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz \
# 	    --o fsaverage_MNI_DKT40+aseg.nii \
# 	    --interp nearest \
# 	    --reg $FREESURFER_HOME/average/mni152.register.dat --inv



################ vcAtlas #################
## vcAtlas is available on the FreeSurfer average surface fsaverage
VC="/projects/crunchie/TOOLBOXES/ATLASES/vcAtlas"
wd="${VC}/NII"
thr=.1

cd $VC
for ll in *.label
do
    echo converting $ll to nii

    $FREESURFER_HOME/bin/mri_label2vol \
	--label $ll \
	--temp $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz \
	--reg $FREESURFER_HOME/average/mni152.register.dat \
	--fillthresh $thr \
	--o $wd/${ll}_thr$thr.nii
done


ribdelta=.1
for ll in *.label
do
    echo converting $ll to nii

    if [[ $ll =~ 'lh' ]]; then
	hem="lh"
    elif [[ $ll =~ 'rh' ]]; then
	hem="rh"
    else
	hem=""
    fi
    echo ....hem=$hem

    $FREESURFER_HOME/bin/mri_label2vol \
	--label $ll \
	--temp $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz \
	--reg $FREESURFER_HOME/average/mni152.register.dat \
	--subject fsaverage \
	--hemi $hem \
	--proj frac 0 1 $ribdelta \
	--fillthresh $thr \
	--o $wd/${ll}_rib${ribdelta}_thr${thr}.nii
done

# $FREESURFER_HOME/bin/mri_annotation2label \
#     --sd $VC/lh.Rosenke_vcAtlas.annot 



################ Place Selectivity #################
CoS="/projects/crunchie/TOOLBOXES/ATLASES/CoS.PlaceSelectivity.Weiner"
wd="${CoS}/NII"
thr=.1

cd $CoS
for ll in *.label
do

    echo converting $ll to nii
    
    $FREESURFER_HOME/bin/mri_label2vol \
	--label $ll \
	--temp $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz \
	--reg $FREESURFER_HOME/average/mni152.register.dat \
	--fillthresh $thr \
	--o $wd/${ll}_thr$thr.nii
done


for ll in *.label
do

    echo converting $ll to nii

    if [[ $ll =~ 'lh' ]]; then
	hem="lh"
    elif [[ $ll =~ 'rh' ]]; then
	hem="rh"
    else
	hem=""
    fi
    echo ....hem=$hem
    
    $FREESURFER_HOME/bin/mri_label2vol \
	--label $ll \
	--temp $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz \
	--reg $FREESURFER_HOME/average/mni152.register.dat \
	--subject fsaverage \
	--hemi $hem \
	--proj frac 0 1 $ribdelta \
	--fillthresh $thr \
	--o $wd/${ll}_rib${ribdelta}_thr${thr}.nii
done

######## DOUBLE-CHECK!! ##########

    freeview -v $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz \
	     $wd/MPM_rh.FG2.label.nii \
	     -l $VC/MPM_rh.FG2.label &


    freeview -v $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz \
	     $wd/MPM_rh.FG2.label_thr.1.nii \
             $wd/MPM_rh.FG2.label_rib.1_thr.1.nii \
	     -l $VC/MPM_rh.FG2.label &


    freeview -v $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz \
	     $wd/MPM_rh.hOc1.label_thr.1.nii \
             $wd/MPM_rh.hOc1.label_rib.1_thr.1.nii \
	     -l $VC/MPM_rh.hOc1.label &
    
##...........

    freeview -v $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz \
	     $wd/rh.CoS.PlaceSelectivity.Weiner.label_thr.1.nii \
             $wd/rh.CoS.PlaceSelectivity.Weiner.label_rib.1_thr.1.nii \
	     -l $CoS/rh.CoS.PlaceSelectivity.Weiner.label &

   mricron $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz &

    ## tkmedit is obsolete -- use freeview
    # To see how well the label is mapped into the functional volume, run




    ## tkmedit is obsolete -- use freeview
    # To see how well the label is mapped into the functional volume, run

    # tkmedit -f $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz \
    # 	    -overlay $wd/MPM_lh.FG1.label.nii \ 
    # 	    -overlay-reg $FREESURFER_HOME/average/mni152.register.dat \
    # 	    -fthresh .5 -fmid 1

# Then load the label with File->Label->LoadLabel. The label should
# overlap with the overlay. The overlap will not be perfect but it
# should be very close.
  
  ## see how good the transform is by loading the orig volume and
  ## applying transform in tkmedit. Then compare the volume with
  ## $subjects_dir/fsaverage/mri/mni305 volume (load it as aux
  ## volume). (Note that you cannot load average_305.mnc which has non
  ## 256^3 volume. tkmedit cannot handle two different sized
  ## volumes. The fsaverage/mri/mni305 is conformed to have 256^3 so
  ## that you can load it under tkmedit.)

  

  
