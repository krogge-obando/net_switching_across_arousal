# 🧠 Network Flexiblity across Arousal State Project

Hello, this repo will store the code needed to reproduce the results in the paper "Arousal state alters network flexibility and moderates its relationship to cognitive task performance"

## ⭐ Highlights

This repo provides code to derive:
-FINDLAB network and parcel time series
-Network flexiblity of the FINDLAB networks and parcels
-code to derive null models used in the paper
-code to compute p-stat from the null models
-generate all the figures from the paper

## 🗒️ How to navigate this repo??

This repo was designed with two folders
1. codes: to see the code used to derive any of the things we stated in our *Highlights* go here, for key information of what the code does go to the section  *codes info*
2. data: the only data provided is the outputs used to derive the figures of this paper. To download the HCP-7T data please go ["here"](https://www.humanconnectome.org/study/hcp-young-adult/data-releases/). For acces to the EEG-fMRI VU data please contact catie.chang@vanderbilt.edu

## Software used in this project

- Matlab version 2024a

Matlab was used to derive the FINDLAB_atlas network and parcel time series, and derive network flexiblity. Please go to external functions to get steps on how to download functions used in *conducting_net_flex.m* and *derive_null_models.m*

- Rstudio version RStudio/2024.12.0+467

Rstudio was used to conduct statistical analysis and visualizations 

## ❗ External functions needed to run the codes

To derive network flexiblity we used the following external functions, please download these codes before attempting to run the code.

1. folder GenLouvain-master at [https://github.com/GenLourvain/HenLouvain](https://github.com/GenLouvain/GenLouvain)
2. flexibilty.m at [http://commdetect.weebly.com/](http://commdetect.weebly.com)

To derive the null models you will need this file

1. brain_benchmark_toolbox-master to get access to this file please email mika.rubinov@vanderbilt.edu

## 💻 Code Info

Due to reduce redundancy we only share the codes needed to run the analysis with the EEG-fMRI VU data.

- **derive_FIND_net_ts.m** - code that conducts dual regression to derive networks time series

- **derive_FIND_parcel_ts.m** - code that derives parcel time series 

- **derive_null_models_ts.m** - code that stores null models

- **run_compute_net_flex.m** code that runs compute_net_flex.m to store the experimental or null models to compute network flexiblity

- **compute_net_flex.m** - function code that derives network flexiblity with either setting for null models or experimental data

- **community_assign_csv.m** - code that stores the community assignment to a csv file for R-studio

- **net_null_model_analysis.R** - code that computes null model analysis to identify p-values

- **parcel_null_model_analysis.R** - code that computes null model analysis to identify p-values

- **violin_plot.R** -code that generates violin plots for this project Figure 2 and 3

- **community_allegiance_dot_plot.R** -code that computes community allegiance and makes the figures for figure 4A

- **community_allegiance_heat_map.R** - code that computes community allegiance as a fraction and makes figures 4C

- **global_flex_task_plots.R** -code that computes and makes figure 5.

## ❓Have Questions

For additional information about the project or how to use the codes feel free to reach out ti my vanderbilt email at kimberly.kundert.obando@vanderbilt.eud
