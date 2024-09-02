#!/bin/bash
#$ -S /bin/bash
#$ -l h_rt=168:00:00,h_vmem=8G,mem_free=8G
#$ -q work.q

# subject ID

CPU=4
ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4

for K in 4031 #4001 4002 4003 4005 4006 4007 4009 4010 4011 4012 4013 4014 4015 4016 4017 4004 4008 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033
do
  ID=${K}

  # set up folders and group space images
  folder=/Volumes/korokdorf/MRPET/coreg_mri/bothsessions/"${ID}"/data/
  # folderCon=/Volumes/korokdorf/MRPET/analysis/bothSessions/1st/"${ID}"/
  folderCon=/Volumes/korokdorf/MRPET/coreg_mri/bothsessions/"${ID}"/data/
  MNI=/Users/alex/Documents/mni_icbm152_t1_tal_nlin_asym_09c.nii
  template=/Volumes/korokdorf/MRPET/mrpet_template.nii.gz

  # -------------- start the transformation -------------- #

  # T1WB -> study-specific template
  # antsRegistration --dimensionality 3 \
  #                  --float 0 \
  #                  --output [${folder}NLantsreg_T1WB_to_template_,${folder}NLantsreg_T1WB_to_template_Warped.nii.gz,${folder}NLantsreg_T1WB_to_template_InverseWarped.nii.gz] \
  #                  --interpolation Linear \
  #                  --winsorize-image-intensities [0.005,0.995] \
  #                  --use-histogram-matching 0 \
  #                  --initial-moving-transform [${template},${folder}T1WB.nii,1] \
  #                  --transform SyN[0.1,3,0] \
  #                  --metric CC[${template},${folder}T1WB.nii,1,4] \
  #                  --convergence [100x70x50x20,1e-6,10] \
  #                  --shrink-factors 4x2x1x1 \
  #                  --smoothing-sigmas 2x1x0.5x0vox \
  #                  --verbose 1

  # meanEPI -> T1WB
  # antsRegistration --dimensionality 3 \
  #                  --float 0 \
  #                  --output [${folder}antscoreg_meanEPI_to_T1WB_,${folder}antscoreg_meanEPI_to_T1WB_Warped.nii.gz,${folder}antscoreg_meanEPI_to_T1WB_InverseWarped.nii.gz] \
  #                  --interpolation Linear \
  #                  --winsorize-image-intensities [0.005,0.995] \
  #                  --use-histogram-matching 0 \
  #                  --initial-moving-transform [${folder}T1WB.nii,${folder}meanEPI.nii,1] \
  #                  --transform Rigid[0.1] \
  #                  --metric MI[${folder}T1WB.nii,${folder}meanEPI.nii,1,32,Regular,0.25] \
  #                  --convergence [1000x500x250x100,1e-6,10] \
  #                  --shrink-factors 4x2x1x1 \
  #                  --smoothing-sigmas 2x1x0.5x0vox \
  #                  --transform Affine[0.1] \
  #                  --metric MI[${folder}T1WB.nii,${folder}meanEPI.nii,1,32,Regular,0.25] \
  #                  --convergence [1000x500x250x100,1e-6,10] \
  #                  --shrink-factors 4x2x1x1 \
  #                  --smoothing-sigmas 2x1x0.5x0vox \

  # T1WB -> meanEPI
  # antsRegistration --dimensionality 3 \
  #                  --float 0 \
  #                  --output [${folder}antscoreg_T1WB_to_meanEPI_,${folder}antscoreg_T1WB_to_meanEPI_Warped.nii.gz,${folder}antscoreg_T1WB_to_meanEPI_InverseWarped.nii.gz] \
  #                  --interpolation Linear \
  #                  --winsorize-image-intensities [0.005,0.995] \
  #                  --use-histogram-matching 0 \
  #                  --initial-moving-transform [${folder}meanEPI.nii,${folder}T1WB.nii,1] \
  #                  --transform Rigid[0.1] \
  #                  --metric MI[${folder}meanEPI.nii,${folder}T1WB.nii,1,32,Regular,0.25] \
  #                  --convergence [1000x500x250x100,1e-6,10] \
  #                  --shrink-factors 4x2x1x1 \
  #                  --smoothing-sigmas 2x1x0.5x0vox \
  #                  --transform Affine[0.1] \
  #                  --metric MI[${folder}meanEPI.nii,${folder}T1WB.nii,1,32,Regular,0.25] \
  #                  --convergence [1000x500x250x100,1e-6,10] \
  #                  --shrink-factors 4x2x1x1 \
  #                  --smoothing-sigmas 2x1x0.5x0vox \
  #                  --verbose 1

  # Apply transformations for T1w to MNI and EPI to MNI
  # antsApplyTransforms -d 3 -v 0 -n BSpline[4] \
  #                     -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_1Warp.nii.gz \
  #                     -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_0GenericAffine.mat \
  #                     -t "${folder}"NLantsreg_T1WB_to_template_1Warp.nii.gz \
  #                     -t "${folder}"NLantsreg_T1WB_to_template_0GenericAffine.mat \
  #                     -i "${folder}"T1WB.nii \
  #                     -r "${MNI}" \
  #                     -o "${folder}"NLantsreg_T1WB_to_MNI.nii

  # antsApplyTransforms -d 3 -v 0 -n BSpline[4] \
  #                     -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_1Warp.nii.gz \
  #                     -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_0GenericAffine.mat \
  #                     -t "${folder}"NLantsreg_T1WB_to_template_1Warp.nii.gz \
  #                     -t "${folder}"NLantsreg_T1WB_to_template_0GenericAffine.mat \
  #                     -t "${folder}"antscoreg_meanEPI_to_T1WB_0GenericAffine.mat \
  #                     -i "${folder}"meanEPI.nii \
  #                     -r "${MNI}" \
  #                     -o "${folder}"NLantsreg_meanEPI_to_MNI_epi2t1.nii.gz

  # antsApplyTransforms -d 3 -v 0 -n BSpline[4] \
  #                     -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_1Warp.nii.gz \
  #                     -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_0GenericAffine.mat \
  #                     -t "${folder}"NLantsreg_T1WB_to_template_1Warp.nii.gz \
  #                     -t "${folder}"NLantsreg_T1WB_to_template_0GenericAffine.mat \
  #                     -t ["${folder}"antscoreg_T1WB_to_meanEPI_0GenericAffine.mat,1] \
  #                     -i "${folder}"meanEPI.nii \
  #                     -r "${MNI}" \
  #                     -o "${folder}"NLantsreg_meanEPI_to_MNI_t12epi.nii.gz

  # Apply transformations for EPI to study-specific template
  # antsApplyTransforms -d 3 -v 0 -n BSpline[4] \
  #                     -t "${folder}"NLantsreg_T1WB_to_template_1Warp.nii.gz \
  #                     -t "${folder}"NLantsreg_T1WB_to_template_0GenericAffine.mat \
  #                     -t "${folder}"antscoreg_meanEPI_to_T1WB_0GenericAffine.mat \
  #                     -i "${folder}"meanEPI.nii \
  #                     -r "${template}" \
  #                     -o "${folder}"NLantsreg_meanEPI_to_template_epi2t1.nii.gz

  # antsApplyTransforms -d 3 -v 0 -n BSpline[4] \
  #                     -t "${folder}"NLantsreg_T1WB_to_template_1Warp.nii.gz \
  #                     -t "${folder}"NLantsreg_T1WB_to_template_0GenericAffine.mat \
  #                     -t ["${folder}"antscoreg_T1WB_to_meanEPI_0GenericAffine.mat,1] \
  #                     -i "${folder}"meanEPI.nii \
  #                     -r "${template}" \
  #                     -o "${folder}"NLantsreg_meanEPI_to_template_t12epi.nii.gz

  # Process contrast images
  if [ "$ID" -eq 4031 ] ; then

    for I in {01..28}
    do
      K=$(printf %04d ${I})
      antsApplyTransforms -d 3 -v 0 -n Linear \
                          -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_1Warp.nii.gz \
                          -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_0GenericAffine.mat \
                          -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz \
                          -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat \
                          -t ["${folder}"antscoreg_T1WB_to_meanEPI_0GenericAffine.mat,1] \
                          -i "${folderCon}"con_${K}_mem.nii \
                          -r "${MNI}" \
                          -o "${folderCon}"con_${K}_mem_mni.nii

      antsApplyTransforms -d 3 -v 0 -n Linear \
                          -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz \
                          -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat \
                          -t ["${folder}"antscoreg_T1WB_to_meanEPI_0GenericAffine.mat,1] \
                          -i "${folderCon}"con_${K}_mem.nii \
                          -r "${template}" \
                          -o "${folderCon}"con_${K}_mem_template.nii
    done

    antsApplyTransforms -d 3 -v 0 -n BSpline[4] \
                      -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_1Warp.nii.gz \
                      -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_0GenericAffine.mat \
                      -t "${folder}"NLreg_T1WB_to_template_1Warp.nii.gz \
                      -t "${folder}"NLreg_T1WB_to_template_0GenericAffine.mat \
                      -t ["${folder}"antscoreg_T1WB_to_meanEPI_0GenericAffine.mat,1] \
                      -i "${folder}"meanEPI.nii \
                      -r "${MNI}" \
                      -o "${folder}"NLantsreg_meanEPI_to_MNI_t12epi_debug.nii.gz

  else

    for I in {01..28}
    do
      K=$(printf %04d ${I})
      antsApplyTransforms -d 3 -v 0 -n Linear \
                          -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_1Warp.nii.gz \
                          -t /Volumes/korokdorf/MRPET/coreg_mri/bothsessions/NLreg_template_to_MNI_0GenericAffine.mat \
                          -t "${folder}"NLantsreg_T1WB_to_template_1Warp.nii.gz \
                          -t "${folder}"NLantsreg_T1WB_to_template_0GenericAffine.mat \
                          -t ["${folder}"antscoreg_T1WB_to_meanEPI_0GenericAffine.mat,1] \
                          -i "${folderCon}"con_${K}_mem.nii \
                          -r "${MNI}" \
                          -o "${folderCon}"con_${K}_mem_mni.nii

      antsApplyTransforms -d 3 -v 0 -n Linear \
                          -t "${folder}"NLantsreg_T1WB_to_template_1Warp.nii.gz \
                          -t "${folder}"NLantsreg_T1WB_to_template_0GenericAffine.mat \
                          -t ["${folder}"antscoreg_T1WB_to_meanEPI_0GenericAffine.mat,1] \
                          -i "${folderCon}"con_${K}_mem.nii \
                          -r "${template}" \
                          -o "${folderCon}"con_${K}_mem_template.nii
    done

  fi

  # ------------------------------------------------------ #

done
