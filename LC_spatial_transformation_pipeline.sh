#!/bin/bash

# subject ID
ID=1001

# set up folders and group space images
folder=/mnt/work/yyi/temp/ED_coreg/"${ID}"/
MNI=/mnt/work/yyi/temp/ED_coreg/mni_icbm152_t1_tal_nlin_asym_09c.nii
template=/mnt/work/yyi/temp/ED_coreg/pilot_template.nii.gz

# LC segmentation
mask=$(ls -t "${folder}"/data/LCmask_"${ID}".nii.gz)

# ---------- prepare images for trasnformation ---------- #

# bias field correct T1 and EPI
N4BiasFieldCorrection -d 3 -v 1 -r 0 -i "${folder}"data/T1mean.nii -o "${folder}"data/T1mean_corrected.nii -s 2 -c [200x150x100x50,1e-6] -b 200
N4BiasFieldCorrection -d 3 -v 1 -r 0 -i "${folder}"data/meanEPI.nii -o "${folder}"data/meanEPI_corrected.nii -s 2 -c [200x150x100x50,1e-6] -b 200

# make an EPI mask (FSL)
/usr/local/fsl/bin/bet "${folder}"data/meanEPI_corrected.nii "${folder}"data/meanEPI_brain -f 0.5 -g 0 -n -m

# resample t1slab to T1 resolution (FreeSurfer)
mri_convert -cs 1 -odt float -rl "${folder}"data/T1mean.nii -rt cubic "${folder}"data/t1slab.nii "${folder}"data/t1slab_1mm.nii

# ------------------------------------------------------ #



# -------------- start the transformation -------------- #

# study template -> MNI
antsRegistrationSyN.sh -d 3 -t s -f "${MNI}" -m "${template}" -o "${folder}"NLreg_template_to_MNI_

# T1 -> study template
antsRegistrationSyN.sh -d 3 -t s -f "${template}" -m "${folder}"data/T1mean_corrected.nii -o "${folder}"data/NLreg_T1mean_to_template_

# T1 -> MNI
antsApplyTransforms -d 3 -v 1 -n BSpline[4] -t "${folder}"NLreg_template_to_MNI_1Warp.nii.gz -t "${folder}"NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"data/NLreg_T1mean_to_template_1Warp.nii.gz -t "${folder}"data/NLreg_T1mean_to_template_0GenericAffine.mat -i "${folder}"data/T1mean_corrected.nii -r "${MNI}" -o "${folder}"data/NLreg_T1mean_to_MNI.nii

# T1 -> EPI
antsRegistrationSyN.sh -d 3 -t r -m "${folder}"data/T1mean_corrected.nii -f "${folder}"data/meanEPI_corrected.nii -x "${folder}"data/meanEPI_brain_mask.nii.gz -o "${folder}"data/coreg_T1mean_to_meanEPI_

# t1slab(resampled) -> T1
antsRegistrationSyN.sh -d 3 -t r -m "${folder}"data/t1slab_1mm.nii -f "${folder}"data/T1mean_corrected.nii -o "${folder}"data/coreg_t1slab_to_T1mean_

# T1 -> MNI
antsApplyTransforms -d 3 -v 1 -n BSpline[4] -t "${folder}"NLreg_template_to_MNI_1Warp.nii.gz -t "${folder}"NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"data/NLreg_T1mean_to_template_1Warp.nii.gz -t "${folder}"data/NLreg_T1mean_to_template_0GenericAffine.mat -i "${folder}"data/T1mean_corrected.nii -r "${MNI}" -o "${folder}"data/NLreg_T1mean_to_MNI.nii

# EPI -> MNI
antsApplyTransforms -d 3 -v 1 -n BSpline[4] -t "${folder}"NLreg_template_to_MNI_1Warp.nii.gz -t "${folder}"NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"data/NLreg_T1mean_to_template_1Warp.nii.gz -t "${folder}"data/NLreg_T1mean_to_template_0GenericAffine.mat -t ["${folder}"data/coreg_T1mean_to_meanEPI_0GenericAffine.mat, 1] -i "${folder}"data/meanEPI_corrected.nii -r "${MNI}" -o "${folder}"data/NLreg_meanEPI_to_MNI.nii

# LC segmentation -> MNI
antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t "${folder}"NLreg_template_to_MNI_1Warp.nii.gz -t "${folder}"/NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"/data/NLreg_T1mean_to_template_1Warp.nii.gz -t "${folder}"/data/NLreg_T1mean_to_template_0GenericAffine.mat -t "${folder}"/data/coreg_t1slab_to_T1mean_0GenericAffine.mat -i "${mask}" -r "${MNI}" -o "${folder}"/data/NLreg_LCmask_to_MNI.nii

# 1st-level stats images -> MNI : ! linear interpolation !
for I in {01..19}
do
	antsApplyTransforms -d 3 -v 1 -n Linear -t "${folder}"NLreg_template_to_MNI_1Warp.nii.gz -t "${folder}"NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"data/NLreg_T1mean_to_template_1Warp.nii.gz -t "${folder}"data/NLreg_T1mean_to_template_0GenericAffine.mat -t ["${folder}"data/coreg_T1mean_to_meanEPI_0GenericAffine.mat, 1] -i "${folder}"data/con_00${I}.nii -r "${MNI}" -o "${folder}"data/con_00${I}_mni.nii
done

# ------------------------------------------------------ #
