maskd="/home/joue/estropharm3/masks_manSeg"

# c3d alEC_PRCpref_left.nii alEC_PRCpref_right.nii pmEC_PHCpref_left.nii pmEC_PHCpref_right.nii -replace 1 2 3 4 -add -o EC_maass.nii

# c3d alEC_PRCpref_left.nii alEC_PRCpref_right.nii -replace 1 2 -add -o alEC_bilat.nii
# c3d pmEC_PHCpref_left.nii pmEC_PHCpref_right.nii -replace 1 2 -add -o pmEC_bilat.nii
# c3d alEC_bilat.nii pmEC_bilat.nii -replace 1 2 -add -o EC_maass.nii
# c3d pmEC_PHCpref_left.nii pmEC_PHCpref_right.nii -replace 3 4 -add -o pmEC_bilat.nii
# c3d alEC_bilat.nii pmEC_bilat.nii -add -o EC_maass.nii
# c3d alEC_bilat.nii pmEC_bilat.nii -replace 1 2 3 4 -add -o EC_maass.nii
# c3d alEC_PRCpref_right.nii pmEC_PHCpref_right.nii -shift 2 -add -o EC_right.nii
# c3d alEC_PRCpref_left.nii pmEC_PHCpref_left.nii -shift 2 -add -o EC_left.nii
# c3d alEC_PRCpref_right.nii pmEC_PHCpref_right.nii -shift 1 -add -o EC_right.nii
# c3d alEC_PRCpref_left.nii pmEC_PHCpref_left.nii -shift 1 -add -o EC_left.nii
# c3d alEC_PRCpref_left.nii pmEC_PHCpref_left.nii -replace 1 2 -add -o EC_left.nii
# c3d alEC_PRCpref_right.nii pmEC_PHCpref_right.nii -replace 1 2 -add -o EC_right.nii

c3d $maskd/al-ErC_L.nii $maskd/pm-ErC_L.nii -replace 1 2 -add -o $maskd/ErC_L.nii
c3d $maskd/al-ErC_R.nii $maskd/pm-ErC_R.nii -replace 1 2 -add -o $maskd/ErC_R.nii
