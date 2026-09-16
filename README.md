# 🧠 Network Switching across Arousal State Project

Hello, this repo will store the code needed to reproduce the results in the paper "Arousal state alters brain network switching and moderates cognitive task performance"

This version of GitHub has been revised 9/16/26 to the current draft of the manuscript that is in review.

## ⭐ Highlights

This repo provides code to derive:
- FINDLAB network and parcel time series
- Network switching of the FINDLAB networks and parcels
- code to compute p-stat from the null models
- generate all the analysis figures from the paper

## 🗒️ How to navigate this repo??

This repo was designed with two folders
1. codes: to see the code used to derive any of the things we stated in our *Highlights* go here, for key information of what the code does go to the section  *codes info*
2. codes/analysis/: to see the codes used to compute the analysis 
3. data: the only data provided is the outputs used to derive the figures 4A and 4C of this paper. To download the HCP-7T data please go ["here"](https://www.humanconnectome.org/study/hcp-young-adult/data-releases/). For acces to the EEG-fMRI VU data please contact catie.chang@vanderbilt.edu

## 🧭 Software used in this project

- Matlab version 2024a

Matlab was used to derive the FINDLAB_atlas network and parcel time series, and derive network flexiblity. Please go to external functions to get steps on how to download functions used in *conducting_net_flex.m* and *derive_null_models.m*

- Rstudio version RStudio/2024.12.0+467

Rstudio was used to conduct statistical analysis and visualizations 

## ❗ External functions needed to run the codes

To derive network switching we used the following external functions, please download these codes before attempting to run the code.

1. folder GenLouvain-master at [https://github.com/GenLourvain/HenLouvain](https://github.com/GenLouvain/GenLouvain)
2. flexibilty.m at [http://commdetect.weebly.com/](http://commdetect.weebly.com)

To get access to the file to derive null models please email mika.rubinov@vanderbilt.edu

## 💻 Code Info

Due to reduce redundancy we only share the codes needed to run the analysis with the EEG-fMRI VU data.

- **derive_FIND_net_ts.m** - code that conducts dual regression to derive networks time series

- **derive_FIND_parcel_ts.m** - code that derives parcel time series 

- **run_compute_net_switch.m** code that runs compute_net_flex.m to store the experimental or null models to compute network flexiblity

- **compute_net_switch.m** - function code that derives network flexiblity with either setting for null models or experimental data

- **community_assign_csv.m** - code that stores the community assignment to a csv file for R-studio

- **net_null_model_analysis.R** - code that computes null model analysis to identify p-values

- **parcel_null_model_analysis.R** - code that computes null model analysis to identify p-values

- **sr_arousal_analysis_w_fig.R** -code that runs all the main analysis for the manuscript and makes the manuscript figures 2 and 3.

- **non_linear_analysis.R** -code that runs the non linear analysis of the manuscript that includes transitions scans and makes figure 4 on the updated manuscript.

- **robust_moderation_tests.R** - code that runs the moderation tests and figure 5.

## ❓Have Questions

For additional information about the project or how to use the codes feel free to reach out to my vanderbilt email at kim.kundert.obando@gmail.com
