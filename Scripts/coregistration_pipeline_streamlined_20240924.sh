#!/bin/bash

path_parent=/mnt/work/yyi/mrpet_both/
MNI=/Users/alex/Documents/mni_icbm152_t1_tal_nlin_asym_09c.nii
template=/Volumes/ALEX3/MRPET/mrpet_template.nii.gz

# first, start off with study-specific template co-registration
antsRegistrationSyN.sh -d 3 -t s -f "${MNI}" -m "${template}" -o "${folder}"NLreg_template_to_MNI_

# subject ID loop starts
for K in 4001 4002 4003 4005 4006 4007 4009 4010 4011 4012 4013 4014 4015 4016 4017 4019 4020 4021 4022 4023 4024 4025 4026 4030 4031 4032 4033 

do

  ID=${K} 


  # set up folders and group space images
  folder=/mnt/work/yyi/mrpet_both/"${ID}"/data/
  folderCon=/Volumes/ALEX3/MRPET/analysis/bothSessions/1st/"${ID}"/ # folder for single volumes of fMRI time-series
  

  # -------------- start the transformation -------------- #

  # T1w -> study-specific template (i.e. zz_template.nii.gz)
  antsRegistrationSyN.sh -d 3 -t s -f "${template}" -m "${folder}"T1WB.nii -o "${folder}"NLreg_T1WB_to_template_

  # T1 slab -> T1WB # only if you have an LC-focused image, such as FLASH or MTw slabs
  antsRegistrationSyN.sh -d 3 -t r -m "${folder}"t1slab.nii -f "${folder}"T1WB.nii -o "${folder}"coreg_t1slab_to_T1WB_

  # T1w -> MNI
  antsApplyTransforms -d 3 -v 0 -n BSpline[4] -t /mnt/work/yyi/mrpet_both/NLreg_template_to_MNI_1Warp.nii.gz -t /mnt/work/yyi/mrpet_both/NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -i "${folder}"T1WB.nii -r "${MNI}" -o "${folder}"NLreg_T1WB_to_MNI.nii



  ########################### !ATTENTION! ##############################
  ############### here, you have two choices. if you have partial-volume fMRIs that has very narrow field of view, try [METHOD A]. otherwise, use [METHOD B] ###############

  # [METHOD A]: 
  # mean functional image (T2*, EPI, representative of your fMRI time series data) -> T1w
  antsRegistrationSyN.sh -d 3 -t r -f "${folder}"meanEPI.nii -m "${folder}"T1WB.nii -o "${folder}"coreg_T1WB_to_meanEPI_

  # mean EPI -> MNI
  antsApplyTransforms -d 3 -v 0 -n BSpline[4] -t "${path_parent}"NLreg_template_to_MNI_1Warp.nii.gz -t "${path_parent}"NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -t ["${folder}"coreg_T1WB_to_meanEPI_0GenericAffine.mat,1] -i "${folder}"meanEPI.nii -r "${MNI}" -o "${folder}"NLreg_meanEPI_to_MNI_epi2t1.nii

  # mean EPI -> study-specific template (uncomment the line below only if you want to do your final analyses in the template space)
  #antsApplyTransforms -d 3 -v 0 -n BSpline[4] -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -t ["${folder}"coreg_T1WB_to_meanEPI_0GenericAffine.mat,1] -i "${folder}"meanEPI.nii -r "${template}" -o "${folder}"NLreg_meanEPI_to_template.nii

  for I in {001..420}
  do
    K=$(printf %04d ${I}) # change this line depending on the numbering format of your single-volume functional image. for example, if your volume is numbered like 'func_000123.nii', it should be $(printf %06d ${I})
    antsApplyTransforms -d 3 -v 0 -n Linear -t "${path_parent}"NLreg_template_to_MNI_1Warp.nii.gz -t "${path_parent}"NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -t ["${folder}"coreg_T1WB_to_meanEPI_0GenericAffine.mat,1] -i "${folderCon}"func_${K}.nii -r "${MNI}" -o "${folderCon}"func_${K}_mni.nii
  done






  # [METHOD B]: 
  # T1w -> mean EPI
  antsRegistrationSyN.sh -d 3 -t r -m "${folder}"meanEPI.nii -f "${folder}"T1WB.nii -o "${folder}"coreg_meanEPI_to_T1WB_

  # mean EPI -> MNI
  antsApplyTransforms -d 3 -v 0 -n BSpline[4] -t "${path_parent}"NLreg_template_to_MNI_1Warp.nii.gz -t "${path_parent}"NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -t "${folder}"coreg_meanEPI_to_T1WB_0GenericAffine.mat -i "${folder}"meanEPI.nii -r "${MNI}" -o "${folder}"NLreg_meanEPI_to_MNI_t12epi.nii

  # mean EPI -> study-specific template (uncomment the line below only if you want to do your final analyses in the template space)
  #antsApplyTransforms -d 3 -v 0 -n BSpline[4] -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -t "${folder}"coreg_meanEPI_to_T1WB_0GenericAffine.mat -i "${folder}"meanEPI.nii -r "${template}" -o "${folder}"NLreg_meanEPI_to_template.nii

  for I in {001..420}
  do
    K=$(printf %04d ${I}) # change this line depending on the numbering format of your single-volume functional image. for example, if your volume is numbered like 'func_000123.nii', it should be $(printf %06d ${I})
    antsApplyTransforms -d 3 -v 0 -n Linear -t "${path_parent}"NLreg_template_to_MNI_1Warp.nii.gz -t "${path_parent}"NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -t "${folder}"coreg_meanEPI_to_T1WB_0GenericAffine.mat -i "${folderCon}"func_${K}.nii -r "${MNI}" -o "${folderCon}"func_${K}_mni.nii
  done



  ########################################################################


  
  for I in {001..420}
  do
    K=$(printf %04d ${I}) # change this line depending on the numbering format of your single-volume functional image. for example, if your volume is numbered like 'func_000123.nii', it should be $(printf %06d ${I})
    antsApplyTransforms -d 3 -v 0 -n Linear -t "${path_parent}"NLreg_template_to_MNI_1Warp.nii.gz -t "${path_parent}"NLreg_template_to_MNI_0GenericAffine.mat -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat -t "${folder}"coreg_meanEPI_to_T1WB_0GenericAffine.mat -i "${folderCon}"func_${K}.nii -r "${MNI}" -o "${folderCon}"func_${K}_mni.nii
  done


  # ------------------------------------------------------ #

done
