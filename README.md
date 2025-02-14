### Brainstem-focused Spatial Transformation Codesets
#### A GitHub repository for the paper: "It is the locus coeruleus! Or… is it? A proposition for analyses and reporting standards for structural and functional magnetic resonance imaging of the noradrenergic locus coeruleus."

---

### Overview
This repository contains the complete set of scripts and methodologies used in our study on advanced MRI techniques focusing on the brainstem nuclei, especially the locus coeruleus. Our goal is to improve the analysis and reporting standards for structural and functional magnetic resonance imaging of this critical area.

### The Paper
- Yeo-Jin Yi, Falk Lüsebrink, Mareike Ludwig, Anne Maaß, Gabriel Ziegler, Renat Yakupov, Michael C. Kreißl, Matthew Betts, Oliver Speck, Emrah Düzel, Dorothea Hämmerer. *Neurobiology of Aging*. 2023.
- DOI: [10.1016/j.neurobiolaging.2023.04.007](https://doi.org/10.1016/j.neurobiolaging.2023.04.007)

### Repository Content
- **Landmark_Insets/** - Anatomical location images delinated on the structural MNI template necessary for inset generation.
- **Scripts/** - Contains all bash and R scripts used for data preprocessing and analysis:
  - **BuildStudySpecificTemplate_AY_20230104.sh** and **create_StudySpecificTemplate_20200201.sh** - Script to build study-specific MRI templates. Uses ANTs.
  - **CalcDistance_FuncLandmarks_1voxPen.m** and **CalcDistance_FuncLandmarks_3voxPen.m** - MATLAB functions for calculating distances in landmark images transformed onto the MNI- or study-specific template space.
  - **LC_spatial_transformation_pipeline.sh** - Main Bash script for executing the spatial transformation pipeline.
  - **coreg_pipeline_antsRegistration.sh** - Bash script for coregistration, using antsRegistration instead of antsRegistrationSyN
  - **coregistration_pipeline_streamlined_20240924.sh** - Streamlined coregistration pipeline script. Uses ANTs.
  - **LandmarkHistogram_singleRater.R** and **LandmarkHistogram_twoRaters.R** - R scripts for generating histogram figures of landmark assessments.
  - **Register_and_Average_LCslabs_20201103.py** - Python script for registering and averaging landmark-based MRI slices.
  - **makeAVIsnapshots.m** - MATLAB script for creating AVI video snapshots from MRI data.
  - **makeHeatmaps.m** - MATLAB script for generating heatmaps from MRI data.
- **Data/** - MRI datasets used in the study (Note: Due to privacy and ethical restrictions, the data is available upon request).
  - **MNI_landmarks_v7_2Labels.nii** - Landmarks drawn on the structural MNI space. This file is required for running the distance calculation MATLAB scripts, **CalcDistance_FuncLandmarks_1voxPen.m** and **CalcDistance_FuncLandmarks_3voxPen.m**.
  - **T1_RegistryCheck_BothApproaches_LQ.mp4** - An example demonstration video of a coregistration method comparison performed on the structural images (T1).
- **Documentation/**
  - **Manual_LCspatialTransformation_20221213.pdf** - Manual describing the procedures and standards proposed in the paper.

### Installation
To successfully run the scripts provided in this repository, specific software installations are required:
- MATLAB (for .m files)
- R (for .R files)
- Python (for .py files)
- Bash (for .sh scripts)
- ANTs, FSL, FreeSurfer (for various image processing scripts)

Below are the quick guides to install each software needed (I use Linus and macOS):

#### MATLAB
- **Step 1**: Visit the [MATLAB](https://www.mathworks.com/products/matlab.html) official website.
- **Step 2**: Choose the appropriate version for your operating system and download it.
- **Step 3**: Follow the installation instructions provided during the setup.

#### R
- **Step 1**: Download R from the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).
- **Step 2**: Select the version for your operating system and follow the provided installation instructions.

#### Python
- **Step 1**: Download Python from the official [Python website](https://www.python.org/downloads/).
- **Step 2**: Select the latest version or a version that your scripts depend on.
- **Step 3**: Install Python by following the setup instructions. Make sure to check 'Add Python to PATH' during installation.

#### Bash
For Windows:
- **Step 1**: Install [Git for Windows](https://gitforwindows.org/), which includes Git BASH.
- **Step 2**: Run the installer and follow the instructions to install Git along with the Bash shell.

For macOS and Linux:
- Bash is typically pre-installed. You can access it through the Terminal application.

#### ANTs (Advanced Normalization Tools)
- **Step 1**: Download Pre-compiled ANTs Binaries. Visit the [ANTs GitHub Releases page](https://github.com/ANTsX/ANTs/releases) to find the latest pre-compiled binaries for your operating system.
- **Step 2**: Download and unzip the binaries to a preferred directory. (My preferred location, as a linux & macOS user, is at /Users/alex/antsbin/.)
- **Step 3**: Add ANTs to PATH:
  ```bash
  export PATH=/path/to/antsbin/bin:$PATH
  ```

#### FSL (FMRIB Software Library)
- **Step 1**: Go to the [FSL website](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation) and choose the appropriate installation instructions for your operating system.
- **tip**: Linux users can often use package managers for easier installation:
  ```bash
  sudo apt-get install fsl
  ```

#### FreeSurfer
- **Step 1**: Register at the [FreeSurfer website](https://surfer.nmr.mgh.harvard.edu/registration.html) to get access to the software.
- **Step 2**: Follow the download and installation instructions provided after registration.
- **Step 3**: Setup environment (add these lines to your shell profile, or simply execute these lines in your terminal window):
  ```bash
  export FREESURFER_HOME=/path/to/freesurfer
  source $FREESURFER_HOME/SetUpFreeSurfer.sh
  ```

### Configuration
After installation, configure the environment paths for Python, ANTs, FSL, and FreeSurfer. This makes sure that these tools are accessible from the command line. Here’s how you can add them to your system’s PATH:

#### For Unix systems (Linux/macOS):
Add the following lines to your `.bashrc` or `.bash_profile`:
```bash
export PATH="/path/to/ANTs:$PATH"
export PATH="/path/to/FSL/bin:$PATH"
export FREESURFER_HOME="/path/to/freesurfer"
source $FREESURFER_HOME/SetUpFreeSurfer.sh
```
#### For Windows:
- **Step 1**: Search for 'Environment Variables' in the Start Menu.
- **Step 2**: Edit the 'Path' variable under System variables to include the paths to Python, ANTs, FSL, and FreeSurfer directories.

Restart your terminal or command prompt to apply these changes. Now, you should be able to run the scripts provided in the repository by navigating to the Scripts directory and executing them as described in the Usage section.

### Usage
To use the scripts:
1. Clone the repository:
   ```
   git clone https://github.com/alex-yi-writes/LC-SpatialTransformation2021.git
   ```
2. Navigate to the script directory:
   ```
   cd LC-SpatialTransformation2021/Scripts
   ```
3. Execute the desired scripts (e.g., for Bash scripts):
   ```
   ./coregistration_pipeline_streamlined_20240924.sh
   ```

### Contributing
I ALWAYS(!) welcome contributions to improve the scripts and methodologies - it's a group effort! :) Please fork the repository and submit a pull request with your proposed changes.

### Contact
For any questions or to request data access, please contact the corresponding author via email at [yeo-jin.yi.15@alumni.ucl.ac.uk](mailto:ucjuyyi@ucl.ac.uk) or [yyi@med.ovgu.de](mailto:yyi@med.ovgu.de).
