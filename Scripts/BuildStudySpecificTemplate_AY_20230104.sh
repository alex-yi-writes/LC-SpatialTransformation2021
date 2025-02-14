#!/bin/bash

# create list of T1w images
for ID in {4001..4030}; 
do
  echo "/mnt/work/yyi/img/$ID/T1w.nii.gz" # adjust paths and filenames here

done > images.txt

# move images.txt to the output directory
mv images.txt /mnt/work/yyi/TEMPLATES/LIDO/output/

# Run from own dir
cd /mnt/work/yyi/TEMPLATES/LIDO/output/

# Generate a group template
antsMultivariateTemplateConstruction2.sh \
  -d 3 -v 6gb -a 1 -i 8 -k 1 -r 1 -n 1 -f 6x4x2x1 \
  -s 4x2x1x0vox -q 100x100x70x20 -t SyN -m CC -c 1 \
  -o zz_ images.txt

# the final output, your study template, will be 'zz_template.nii.gz'
