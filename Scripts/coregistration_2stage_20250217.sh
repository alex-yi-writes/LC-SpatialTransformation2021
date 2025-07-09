#!/bin/bash
# a small disclaimer from yeo-jin alex:
# 	this pipeline performs only the nonlinear anatomical-functional transformation part! please try this script when there are persistently difficult transformation that ever so slightly doesn't line up. heads up: this pipeline takes a lot of resources and time. but sometimes you just gotta do what you gotta do...

for ID in sub-201 sub-202
do

	# ===== set paths and vars ===== #

	# paths
	path_parent=/Users/alex/paperwriting/1315/data/segmentation/
	path_anat=/Users/alex/paperwriting/1315/data/segmentation/${ID}v1s1/anat
	path_func=/Users/alex/paperwriting/1315/data/segmentation/${ID}v1s1/func
	path_trx=/Users/alex/paperwriting/1315/data/segmentation/${ID}v1s1/transx
	path_roi=/Users/alex/paperwriting/1315/data/segmentation/${ID}v1s1/anat/roi # in this folder, you should already have two binary ROI masks, e.g. brainstem or mPFC mask in my case, that each in functional and anatomical space. these will be used in the second-stage nonlinear transformation. you can use any of your favourite segmentation tools. i used fastsurfer (https://deep-mi.org/research/fastsurfer/) for creating them.

	# anatomical images
	t1w=${path_anat}/t1/m${ID}v1s1_run-01_T1w_0pt35 # t1w, resampled

	# mean functional images
	func_origenc=${path_func}/meanu${ID}v2s2_task-origenc_run-03_bold

	# prepare workspace
	mkdir ${path_roi}
	mkdir ${path_trx}


	# ======== strip skulls ======== #
	# mri_synthstrip is a part of freesurfer, but it exists as a standalone suit here: https://surfer.nmr.mgh.harvard.edu/docs/synthstrip/

	mri_synthstrip -i ${t1w}.nii.gz -o ${t1w}_stripped.nii.gz
	mri_synthstrip -i ${func_origenc}.nii.gz -o ${func_origenc}_stripped.nii.gz

	#======== nonlinear anat-func transformation ======== #

	# -- whole-brain to functional image

	# 1) i used rigid+affine+syn because distortion correction of my functional image was not optimal. if you have a functional image that has very little distortion and aligns well with the corresponding anatomical image, please use rigid only to save time.
	antsRegistrationSyN.sh -d 3 -n 4 -t s -m ${t1w}_stripped.nii.gz -f ${func_origenc}_stripped.nii.gz -o ${path_trx}/wb2origenc_

	# 2) now, we do the second stage of the transformation.
	# hint: you can also try another cost function : -m "CC[${func_origenc}_stripped.nii.gz,${t1w}_stripped.nii.gz,1,3]"
	antsRegistration -d 3 \
	  -o "[${path_trx}/BrainStemref_origenc_,${path_trx}/BrainStemref_origenc_Warped.nii.gz]" \
	  --initial-moving-transform ${path_trx}/wb2origenc_0GenericAffine.mat \
	  --initial-moving-transform ${path_trx}/wb2origenc_1Warp.nii.gz \
	  -t "SyN[0.015,3,0]" \
	  -m "MI[${func_origenc}_stripped.nii.gz,${t1w}_stripped.nii.gz,0.7,64,Regular,0.5]" \
	  -c "[1000x500x250x100,1e-8,15]" \
	  -s 3x1.5x1x0 \
	  -f 6x3x2x1 \
	  -x "[${path_roi}/regmask/BrainStem_on_func_origenc.nii.gz,${path_roi}/regmask/BrainStem_on_anat.nii.gz]"

	# ~~~ now, you can incorporate the new warps to the previous pipeline and see if it worked. always check the outputs of previous steps first before applying the transformations to the contrast images, because this method doesn't always work.
	antsApplyTransforms -d 3 -v 1 -n Linear \
		-r /Users/alex/templates/mni_icbm152_t1_tal_nlin_asym_09c.nii.gz \
		-i ${func_origenc}.nii.gz \
		-o ${path_trx}/func_origenc_on_mni.nii.gz \
		-t ${path_parent}/NLreg_template_to_MNI_1Warp.nii.gz \
		-t ${path_parent}/NLreg_template_to_MNI_0GenericAffine.mat \
		-t ${path_trx}/NLreg_T1WB_to_template_1Warp.nii.gz \
		-t ${path_trx}/NLreg_T1WB_to_template_0GenericAffine.mat \
		-t ${path_trx}/BrainStemref_origenc_2InverseWarp.nii.gz \
  		-t ${path_trx}/wb2origenc_1InverseWarp.nii.gz \
		-t [${path_trx}/wb2origenc_0GenericAffine.mat,1] \

done