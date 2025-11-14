# Inter-individual similarities in internal models of the world shape similarities in the perception and neural processing of scenes.

This repository contains the code and data for the behavioural component of the paper:
**Engeser & Kaiser (2025). Inter-individual similarities in internal models of the world shape similarities in the perception and neural processing of scenes.**

It includes experimental code, behavioural data analysis scripts, and materials for studying how individual differences in world models shape scene perception.

## Data
Large model files (e.g., VGG16 network and NIfTI files containing fMRI data) are available on [OSF](https://osf.io/zjtwx/).
Download and place them in the appropriate folder.

Basic structure 
pep_wp4_beh-fMRI/

│

├── fMRI

│ ├── code

│ ├── derivatives

│ ├── drawings

│ ├── localizer

│ ├── MNI_ROIs

│ ├── photos

│ ├── sourcedata

│ ├── stimuli

│ ├── vgg16_imagenet

│ └── vgg16_places265

│

├── behavior

│ ├── experiment/ # Behavioural experiments (Matlab + Psychtoolbox)

│ │ ├── data

│ │ ├── functions

│ │ ├── instructions

│ │ ├── stimuli

│ │ ├── trial_matrices

│ │ └── Run_wp4_beh.m

│ │ 

│ ├── analysis/ # Data analysis scripts (Matlab, R, Python)

│ │ ├── functions

│ │ ├── vgg16_imagenet

│ │ ├── vgg16_places265

│ │ ├── analysis_script_exp1.mlx

│ │ └── analysis_script_exp2.mlx

│ │ 

│ ├──image_similarities/ # Drawings and generated photorealistic images

│ │ ├── drawings_draw3D

│ │ ├── drawings_human_rated

│ │ └── pictures

│

└── README.md



[...]

Under experiment, you can find the code and stimuli for the behavioural experiments, including the categorization tasks (Experiments 1 and 2) and the rating task (Experiment 1). The experiments are implemented in MATLAB using Psychtoolbox-3.

Under analysis, you can find the code for the data analysis. Each experiment has a dedicated live script that controls the analysis and calls the respective functions from the functions folder.
For the analysis of image similarity using deep neural networks (DNNs), please download the required DNN files from OSF (https://osf.io/zjtwx/) and place them in the appropriate folder. The analysis is primarily conducted in MATLAB, with R used for linear mixed-effects (LME) modelling and Python for CLIP and DINO feature extraction.

Under image_similarities, you can find the participant drawings and the corresponding photorealistic images created from these drawings.
This folder also contains the code for the behavioural drawing similarity rating experiment, implemented in PsychoPy.

## Citation
If you use this repository, please cite:

Engeser, L., & Kaiser, D. (2025). *Inter-individual similarities in internal models of the world shape similarities in the perception and neural processing of scenes.*

## Contact
For questions or collaboration:
- **Lead author:**  Micha Engeser
- **Email:** michaengeser[at]gmail.com







Code and Data will be added here in a user-friendly repository combining all code of the study

## In the meantime check out repos of project parts
https://github.com/michaengeser/pep_wp4 - for behavioural study

https://github.com/michaengeser/pep_wp4_fMRI - for fMRI study


## Contact
For questions or collaboration:
- **Lead author:**  Micha Engeser
- **Email:** michaengeser[at]gmail.com
