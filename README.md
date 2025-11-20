# Inter-individual similarities in internal models of the world shape similarities in the perception and neural processing of scenes.

This repository contains the code and data for the behavioural component of the paper:
**Engeser & Kaiser (2025). Inter-individual similarities in internal models of the world shape similarities in the perception and neural processing of scenes.**

It includes experimental code, behavioural data analysis scripts, and materials for studying how individual differences in world models shape scene perception.

## Additional Data Links
Download and place the data files in the appropriate folder (see below for folder structure).

VGG16 network files are available on [OSF](https://osf.io/zjtwx/).

Functional and anatomical sourcedata files can be found here: >data not uploaded yet< (no raw anatomical scans are shared to protect participants' privacy)

Preprocessed files of task runs are stored here: >data not uploaded yet< (unzip and put in fMRI/derivatives folder)

Preprocessed files of localizer runs are stored here: >data not uploaded yet< (unzip and put in fMRI/derivatives folder)



## Basic structure 
pep_wp4_beh-fMRI/

│

├── **fMRI** (contains code for fMRI experiment and analysis)

│ ├── code

│ │ ├── analysisPipeline.mlx (main analysis script for fMRI experiment)

│ │ ├── run_experiment_fMRI.m (experiment script)

│ │ └── utilities

│ │ 

│ ├── derivatives (unzip and place downloaded preprocessed files here!)

│ ├── drawings

│ ├── localizer

│ ├── MNI_ROIs

│ ├── photos (private photos are not shared, but IS-RDM is provided to replicate analysis)

│ ├── sourcedata (unzip and place downloaded raw files here!)

│ ├── stimuli

│ ├── vgg16_imagenet (put vgg16_imagenet model from OSF here)

│ └── vgg16_places265 (put vgg16_places265 model from OSF here)

│

├── **behavior** (contains code for behavioral experiment and analysis)

│ ├── experiment

│ │ ├── data

│ │ ├── functions

│ │ ├── instructions

│ │ ├── stimuli

│ │ ├── trial_matrices

│ │ └── Run_wp4_beh.m

│ │ 

│ ├── analysis

│ │ ├── functions

│ │ ├── vgg16_imagenet (put vgg16_imagenet model from OSF here)

│ │ ├── vgg16_places265 (put vgg16_places265 model from OSF here)

│ │ ├── analysis_script_exp1.mlx (main analysis script for first behavioral experiment)

│ │ └── analysis_script_exp2.mlx (main analysis script for second behavioral experiment)

│ │ 

│ ├──image_similarities

│ │ ├── drawings_draw3D (generated images)

│ │ ├── drawings_human_rated (raw drawings and script for human rating experiment)

│ │ └── pictures (private photos are not shared, but IS-RDM is provided to replicate analysis)

│

└── **README.md**

The scripts to run the experiments are implemented in MATLAB and Psychtoolbox-3.
Most analysis is also coded in MATLAB, but R is used for linear mixed-effects (LME) modelling and Python for CLIP and DINO feature extraction.

To reproduce the analysis, download the respective data (see above) and run the MATLAB live scripts (.mlx). These are markdown scripts that should be self-explanatory. If anything is unclear, please contact me.

## Citation
If you use this repository, please cite:

Engeser, L., & Kaiser, D. (2025). *Inter-individual similarities in internal models of the world shape similarities in the perception and neural processing of scenes.*

## Contact
For questions or collaboration:
- **Lead author:**  Micha Engeser
- **Email:** michaengeser[at]gmail.com


