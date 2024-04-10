## Overview
This repository contains two folders, **scripts** and **models**.

**models** contains MATLAB code implementing Models 1-5 described in the manuscript under Supplemental Information - Description of calcium binding models. This folder also contains a subfolder **model_fits** which includes fitted parameter values for each model as well as the dataset used to fit the model.

**sciprts** contains two scripts. ``plot_all_results`` reproduces all the quantitative results presented in the paper and ``example_fitting_script`` provides example code for fitting a model to data.


## System requirements

Running the scripts requires only a standard computer. This software has been tested using MATLAB 2023a running on Windows 10 and should be compatbile with Linux and MacOS.

Required MATLAB toolboxes
+ Statistics and Machine Learning Toolbox
+ Optimization Toolbox
+ Global Optimization Toolbox

## Installation guide
```
git clone https://github.com/niklasbrake/AMPAR_permeation_modelling
cd AMPAR_permeation_modelling/scripts
matlab
```
This should run in less than a minute and open a matlab instance with the working directory set to the **scripts** folder.

## Demo
An example of the model fitting procedure can be run with the script
```
example_fitting_script
```
This script performs bootstrap samples of the calcium block data and refits model 4.2 to assess uncertainty in the model's parameter values. Running this script on a standard laptop takes approximately 20 min for 100 bootstrap sample. 

After the script runs, the matlab struct variable, ``P_bootstrap``, will contain four fields, **A1_A2**, **A1_A2G2**, **A1_A2G2_C3**, and **A1_A2_R607EG2** corresponding to each of the four AMPARs, each containing an Nx9 array of parameter values, where N corresponds to the number of bootstrap sampling. The nine parameters are in the following order, k1*, k1r*, k2*, k2r*, k3*, k4r*, k4*, d1, d3.


## Instructions for use
Quantiative results in our manuscript can be reproduced by running the script
```
plot_all_results
```
This script should produce 12 plots that reproduce the modelling results in Extended Data Figs. 8-10.

## License
<a rel="license" href="https://www.gnu.org/licenses/gpl-3.0.html#license-text"><img alt="GPL3" style="border-width:0" src="https://www.gnu.org/graphics/gplv3-127x51.png" /></a><br />This repository is licensed under a <a rel="license" href="https://www.gnu.org/licenses/gpl-3.0.html#license-text">GNU GENERAL PUBLIC LICENSE Version 3</a>.
