# Cross-Lagged Model

## Description

## Table of Contents
- [Getting Started](#getting-started)
- [Usage](#usage)
- [File Descriptions](#file-descriptions)

## Getting Started

1. Clone the repository.
2. Install the required dependencies for running the preprocessing. These dependencies can be found in `config.txt`.
3. Navigate to `preprocess.py` for running the preprocessing.
4. Nagivate to `networkanalysis.R` for running the Cross-lagged Model. The R dependencies will be automaticaly installed when running the network object.
4. Read the documentation.
5. Change the config if nessecary.

After running the `networkanalysis.R`, multiple plots will be visable in the given image directory.
The analysis can take some time depening om the hardware and the amount of selected boots.
When changing the plot parameters it's best to build and set the coefficient manualy, so that the calculation only has to be done 1 time.

## Usage

`networkanalysis.R` is the main file to use to calculate the cross-lagged model. This will be done using the glmnet package as a way to calculate the edge strength for the different timesteps.
the function `plot_functions()` will calculate the coefficient for you and uses that to make all the plots.

For manual use:
1. Call the `build_coef()` function.
2. Call the plot functions you wan to use.
If you want to change the plot parameters I would advise to first call the `build_coef()` once and save the output to a variable. With the `set_coef()` function you can save this coefficient everytime you want to make changes to the code. This way it will save a lot of time.

## File Descriptions
- `src/preprocess.ipynb`: Preprocesses the data, such as removing NaN values and enabeling the covariances.
- `src/calculate_descriptive.ipynb`:  Different way to calculate the descriptives from the dataset. This is done because the function to calculate the descriptive in R uses the data imputed by the covriance.
- `src/networkanalysis.R`: Main file for calculating the Cross-Lagged Model, using the glmnet package.

## Author
Code is writen by Tycho Stam
https://github.com/kingilsildor
