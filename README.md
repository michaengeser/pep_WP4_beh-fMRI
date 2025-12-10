# Inter-individual similarities in internal models of the world shape similarities in the perception and neural processing of scenes

This repository contains the code and data of the paper:

### Engeser & Kaiser (2025). Inter-individual similarities in internal models of the world shape similarities in the perception and neural processing of scenes. (Preprint: https://doi.org/10.31234/osf.io/hd6je_v1)

It includes experimental code, behavioral analysis scripts, and materials for investigating how individual differences in world models influence scene perception.


## Additional Data Links

Some large files are hosted externally and must be downloaded separately.

VGG16 network files: available on OSF → https://osf.io/zjtwx/

Functional and anatomical sourcedata files: https://zenodo.org/records/17653128
(Raw anatomical scans are not shared to protect participant privacy.)

Preprocessed task-run files: https://zenodo.org/records/17661600
(Unzip and place in the fMRI/derivatives folder.)

Preprocessed localizer-run files: https://zenodo.org/records/17701429
(Unzip and place in the fMRI/derivatives folder.)


To reproduce the figures from the paper, the intermediate files included in this repository are sufficient.
Download the fMRI data only if you wish to reproduce the full preprocessing or time-course extraction. Corresponding code is included in the repository.


## Repository Structure

pep_wp4_beh-fMRI/

│

├── fMRI/                     # fMRI experiment + analysis code

│   ├── code/

│   │   ├── analysisPipeline.mlx     # Main analysis script (fMRI)

│   │   ├── run_experiment_fMRI.m    # fMRI experiment script

│   │   └── utilities/

│   │

│   ├── derivatives/          # Place preprocessed files here (after unzipping)

│   ├── drawings/

│   ├── localizer/

│   ├── MNI_ROIs/

│   ├── photos/               # Private photos not shared; IS-RDM included

│   ├── sourcedata/           # Place raw data here (after unzipping)

│   ├── stimuli/

│   ├── vgg16_imagenet/       # Add VGG16 ImageNet model from OSF

│   └── vgg16_places265/      # Add VGG16 Places model from OSF

│

├── behavior/                 # Behavioral experiment + analysis code

│   ├── experiment/

│   │   ├── data/

│   │   ├── functions/

│   │   ├── instructions/

│   │   ├── stimuli/

│   │   ├── trial_matrices/

│   │   └── Run_wp4_beh.m

│   │

│   ├── analysis/

│   │   ├── functions/

│   │   ├── vgg16_imagenet/       # Add VGG16 ImageNet model from OSF

│   │   ├── vgg16_places265/      # Add VGG16 Places model from OSF

│   │   ├── analysis_script_exp1.mlx   # Analysis for Experiment 1

│   │   └── analysis_script_exp2.mlx   # Analysis for Experiment 2

│   │

│   └── image_similarities/

│       ├── drawings_draw3D/        # Generated images

│       ├── drawings_human_rated/   # Raw drawings + human rating script

│       └── pictures/               # Private images not shared; IS-RDM included

│

└── README.md


## How to Use This Repository
### Experiments 
(The behavioral and fMRI experiments are implemented in MATLAB with Psychtoolbox-3). 

Run_wp4_beh.m for the behavioral experiment

run_experiment_fMRI.m for the fMRI experiment




### Analysis
(Most analyses are implemented in MATLAB, with R for linear mixed-effects (LME) modeling Python for CLIP and DINO feature extraction).

To reproduce the analyses:

Download the required data (see links above) and place files in the indicated folders.

Run the MATLAB live scripts (.mlx), which include step-by-step explanations.

If anything is unclear, feel free to contact me.



## Citation

If you use this repository, please cite:

Engeser, M., & Kaiser, D. (2025). Inter-individual similarities in internal models of the world shape similarities in the perception and neural processing of scenes.



## Contact

For questions or collaborations:

Lead author: Micha Engeser

Email: michaengeser[at]gmail.com
