import os
from nipype.interfaces.ants import RegistrationSynQuick, AverageImages

subject_ids = range(4001, 4034)
base_path = "/mnt/work/yyi/LCslabs/"

for subject_id in subject_ids:
    subject_path = os.path.join(base_path, str(subject_id))
    ref_image = os.path.join(subject_path, "MTC1.nii.gz")  # ref image is always the first acquisition just because

    if not os.path.exists(ref_image):
        print(f"Reference image MTC1 missing for subject {subject_id}")
        continue

    registered_files = [ref_image] 

    # Perform registration
    for i in range(2, 5):
        moving_image_path = os.path.join(subject_path, f"MTC{i}.nii.gz")
        if os.path.exists(moving_image_path):
            reg = RegistrationSynQuick()
            reg.inputs.fixed_image = ref_image
            reg.inputs.moving_image = moving_image_path
            reg.inputs.output_prefix = os.path.join(subject_path, f"regMTC{i}_")
            reg.inputs.args = '-t r'  # always rigid, don't change
            reg_out = reg.run()
            registered_files.append(reg_out.outputs.warped_image)
        else:
            print(f"moving image MTC{i} missing for subject {subject_id}")

    # Average the registered images
    if len(registered_files) > 1:
        avg = AverageImages()
        avg.inputs.dimension = 3
        avg.inputs.output_average_image = os.path.join(subject_path, "LCslab_averaged.nii.gz")
        avg.inputs.normalize = False
        avg.inputs.images = registered_files
        avg.run()
    else:
        print(f"not enough images to average for subject {subject_id}")
