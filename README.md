# Inter-individual similarities in internal models predict similarities in the perception and neural processing of scenes

This repository contains the code and data of the paper:

### Engeser & Kaiser (2026). Inter-individual similarities in internal models predict similarities in the perception and neural processing of scenes. (Preprint: https://osf.io/preprints/psyarxiv/hd6je_v3)

It includes the experimental code as well as the analysis scripts and materials for the inter-subject representational similarity analysis.

## Installation 
Download or clone the repository (~7GB without fMRI data). Download VGG16 network files (https://osf.io/zjtwx/) and place them in the folders named accordingly. 

To reproduce the outputs from the paper, the intermediate files included in this repository are sufficient.
Download the fMRI data only if you wish to reproduce the full preprocessing or time-course extraction. Corresponding code is included in the repository.

Download and installation should take <30 min. on a "normal" computer (without fMRI data).

The code is mostly written in MATLAB (R2022a) on Windows 11. Installation on other operating systems might need minor adjustments.

Additionally, install the following software dependencies: Psychtoolbox (3.0.19) for running the experiment and SPM12 for fMRI analysis

## Demo
As a demo, you can run the analysis live script of the second experiment (pep_WP4_beh-fMRI\behavior\analysis\analysis_script_exp2.mlx). This contains the core steps of the inter-subject RSA framework, including building IS-RDMs of the drawings and task performance, comparing the IS-RDMs, and statistical testing. More detailed information on running the code can be found in the Markdown sections of the live script.

The demo should take <10 min. on a "normal" computer.

## Data links
Some large files are hosted externally and must be downloaded separately.

VGG16 network files: available on OSF → https://osf.io/zjtwx/

Functional and anatomical sourcedata files: https://zenodo.org/records/17653128
(Unzip and place in the fMRI/sourcedata folder. Raw anatomical scans are not shared to protect participant privacy.)

Preprocessed task-run files: https://zenodo.org/records/17661600
(Unzip and place in the fMRI/derivatives folder.)

Preprocessed localizer-run files: https://zenodo.org/records/17929617
(Unzip and place in the fMRI/derivatives folder.)


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


### Experiments 
The behavioral and fMRI experiments are implemented in MATLAB (R2022a - 9.12.0.2327980) with Psychtoolbox (3.0.19).

Run_wp4_beh.m for the behavioral experiment

run_experiment_fMRI.m for the fMRI experiment


### Analysis
Most analyses are implemented in MATLAB (R2022a - 9.12.0.2327980), including SPM12. Additionally, R (4.3.1) was used for linear mixed-effects (LME) modeling, and Python (3.10) for accessing CLIP, DINO, Gemini, and MPNet.

To reproduce the analyses:

Download the required data and material (see links above) and place files in the indicated folders.

Run the MATLAB live scripts (.mlx), which include step-by-step explanations.

If anything is unclear, feel free to contact me.



## Citation

If you use this repository, please cite:

Engeser, M., & Kaiser, D. (2026). Inter-individual similarities in internal models predict similarities in the perception and neural processing of scenes.



## Contact

For questions or collaborations:

Lead author: Micha Engeser

Email: michaengeser[at]gmail.com
