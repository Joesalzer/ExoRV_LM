This repository includes code, data, and summaries corresponding to the paper "Searching for Low-Mass Exoplanets Amid Stellar Variability with a Fixed Effects Linear Model of Line-by-Line Shape Changes".

The code uses R (https://www.r-project.org/) and assumes you have installed the following packages:
tidyverse, rhdf5, Matrix, patchwork, collapse, parallel, pbmcapply

Data used for this paper can be found at: ZENODO LINK. 

If you come across any issues or bugs, please contact Joseph Salzer at jsalzer@wisc.edu.

Below are a list of all scripts in this repo with a short explanation of what they do. There is also a minimal working example to re-produce one of the models that we fit.

file: readSpectra.R <br>
> r script that contains helper functions for use throughout this project

file: rv_lm.R <br>
> Fit a two-way fixed effects linear model, with options to include covariates and which csv file to be used.

> This script has some variables that need to be set by the user (such as the working directory, what the response column is called, what the time column is called, etc.) inside a block with the heading "## variables to be set by the user ##". This script assumes that the covariates are Gaussian fit parameters like depth, width, a and b (these begin with "fit_gauss") or the projected HG coefficients (these begin with "proj_hg_coeff_"). These are only use for naming the model, whereas any custom covariates can be used but will result in the model name being "Gauss=none_HG=none". Feel free to edit this naming convention using lines 55-64 under "# name the current model, using Gaussian fit parameters and hg covariateNames".

> To run the model, go into the terminal at the working directory and run: Rscript rv_lm.R "CSV FILE" "VAR1,VAR2,..." If everything runs smoothly a folder called "models" will be created in your working directory with a subfolder of the model's name and metrics that have to do with the fit model itself as a file called "model.rds".  If you want to run the model without any covariates, leave the variables blank. 

file: rv_lm_leverages.R <br>
> This file calculates the leverages for the design matrix of one of the models. This is done in batches to speed up computation. To run: Rscript rv_lm_leverages.R "MODEL NAME" "BATCH SIZES". 

> This script has some variables that need to be set by the user (such as the working directory, what the response column is called, what the time column is called, etc.) inside a block with the heading "## variables to be set by the user ##". 

file: rv_boot.R <br>
> This script creates the wild bootstrap samples for the uncertainty quantification of the cleaned RVs using a model fit by rv_lm.R. To run: Rscript rv_boot.R "MODEL NAME" "NUMBER OF BOOTS". 

> This script has some variables that need to be set by the user (such as the working directory, what the response column is called, what the time column is called, etc.) inside a block with the heading "## variables to be set by the user ##". 

file: rv_cv.R <br>
> This script runs the cross-validation (LOODCV, LOOWCV, LOOMCV) procedure of a given model as outlined by the paper. To run: Rscript rv_cv.R "MODEL NAME" "BLOCK SIZE". Note that the CV procedures for the LASSO model were fit using the modeling.Rmd markdown script, but the underlying code is the same as this script.

> This script has some variables that need to be set by the user (such as the working directory, what the response column is called, what the time column is called, etc.) inside a block with the heading "## variables to be set by the user ##". 

file: clean_data.Rmd <br>
> This markdown file uses the directory "line_property_files" to join all of the absorption line data together. The output is completeLines.csv, which is a data file for every line's shape properties, line-specific information, and contaminated radial velocity. This script also provides some summary information and checks of individual lines (like lines that have constant b1 and width). Some lines and days are removed from the analysis (and some lines failed to be read-in) as specified in the paper. Projected HG coefficients are also calculated, reading in the file project_template/project_template_deriv_onto_gh.h5
	
file: eda.Rmd <br>
> This markdown file contains the EDA for this project, including visuals before models were fit. This file creates some of the figures from the paper: the heatmaps (figure 1) and the RV/depth of select lines (figure 2)

file: modeling.Rmd <br>
> This markdown file  is a testing grounds for all of our modeling, including the full-data analysis, bootstrap, and cv. This script includes the fitting of the Common Slopes and Full Model w/ LASSO.

file: results.Rmd <br>
> This markdown file provides key results from the full-data analysis, bootstrap, and cv from the models produced by rv_lm.R/rv_boot.R/rv_cv.R. This script also creates visuals for after the modeling was done: figure 3, 4 and 5.
	
file: pcaRV.Rmd <br>
> This markdown file  does PCA on the RV time series of lines. This was a preliminary analysis that revealed that around 250 of the analyzed absorption lines had a bias before and after the fire.


*Minimal working example*

This minimal working example will reproduce the results from the Full Model of the paper.

1. Make sure that you have the following packages downloaded: Matrix, tidyverse, stringr, rhdf5.
2. Download the following files: completeLines.csv, rv_lm.R, readSpectra.R.
3. Put all files into the same working directory. In the file "rv_lm.R" replace the WD_DATA object with this working directory and save the file.
4. Go to terminal at this directory, and run: Rscript rv_lm.R completeLines.csv fit_gauss_a,fit_gauss_b,fit_gauss_depth,fit_gauss_sigmasq,proj_hg_coeff_0,proj_hg_coeff_2,proj_hg_coeff_3,proj_hg_coeff_4,proj_hg_coeff_5,proj_hg_coeff_6,proj_hg_coeff_7,proj_hg_coeff_8,proj_hg_coeff_9,proj_hg_coeff_10
5. If everything runs smoothly a folder called "models" will be created in your working directory with a subfolder of the model's name and metrics that have to do with the fit model itself as a file called "model.rds"
6. Results can be viewed by opening this model in R via the readRDS function.