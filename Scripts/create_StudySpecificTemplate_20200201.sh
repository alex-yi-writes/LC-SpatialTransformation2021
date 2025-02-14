#!/bin/bash

CPUs=4

# write a list of WBT1
find . -name '*nii' -o -name '*.nii' > mprage.txt

# generate a group template
antsMultivariateTemplateConstruction2.sh -d 3 -a 1 -i 4 -k 1 -r 1 -f 6x4x2x1 -s 4x2x1x0vox -q 100x100x70x20 -t SyN -m CC -c 0 -o mprage.txt
