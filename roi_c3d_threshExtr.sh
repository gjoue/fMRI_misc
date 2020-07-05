  

  # volume stats
  c3d e3_ErC.nii.gz -lstat
  c3d mean_wmT1.nii e3_ErC.nii.gz -lstat
#  c3d e3_ErC.nii.gz -split
  c3d mean_wmT1.nii e3_seg_ErC.nii.gz -lstat >> e3_seg_ErC_VOLS.txt

  c3d e3_seg_ErC.nii.gz -split -oo al-ErC_L.nii pm-ErC_L.nii al-ErC_R.nii pm-ErC_R.nii

  c3d e3_seg_ErC.nii.gz -thresh 0.4 1.4 1 0 -o al-ErC_L.nii
  c3d e3_seg_ErC.nii.gz -thresh 1.4 2.4 1 0 -o pm-ErC_L.nii
  c3d e3_seg_ErC.nii.gz -thresh 2.4 3.4 1 0 -o al-ErC_R.nii
  c3d e3_seg_ErC.nii.gz -thresh 3.4 4.4 1 0 -o pm-ErC_R.nii

  c3d al-ErC_L.nii -resample-mm 2.0x2.0x2.0mm -o al-ErC_L.2mm.nii
  c3d pm-ErC_L.nii -resample-mm 2.0x2.0x2.0mm -o pm-ErC_L.2mm.nii
  c3d pm-ErC_R.nii -resample-mm 2.0x2.0x2.0mm -o pm-ErC_R.2mm.nii
  c3d al-ErC_R.nii -resample-mm 2.0x2.0x2.0mm -o al-ErC_R.2mm.nii
