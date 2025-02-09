file: readSpectra.R <br>
	 -r script that contains helper functions for use throughout this project

file: clean_data.Rmd <br>
	-Uses the directory "line_property_files" to join all of the absorption line data together. The output is completeLines.csv, which is a data file for every line's shape properties, line-specific information, and contaminated radial velocity. This script also provides some summary information and checks of individual lines (like lines that have constant b1 and width). Some lines and days are removed from the analysis (and some lines failed to be read-in) as specified in the paper. Projected HG coefficients are also calculated, reading in the file project_template/project_template_deriv_onto_gh.h5
	
file: eda.Rmd <br>
    EDA for this project, includes visuals before models were fit. This file creates some of the figures from the paper: the heatmaps (figure 1) and the RV/depth of select lines (figure 2)

file: pcaRV.Rmd <br>
    This file does PCA on the RV time series of lines. This was a preliminary analysis that revealed that around 250 of the analyzed absorption lines had a bias before and after the fire.

file: modeling.Rmd <br>
    Testing grounds for all of our modeling, including the full-data analysis, bootstrap, and cv. This script includes the fitting of the Common Slopes and Full Model w/ LASSO.

file: results.Rmd <br>
    Provides some key results from the full-data analysis, bootstrap, and cv along with visuals for after the modeling was done
		
file: rv_lm.R <br>
> Fit a two-way fixed effects linear model, with options to include covariates and which csv file to be used.

> This script has some variables that need to be set by the user (such as the working directory, what the response column is called, what the time column is called, etc.) inside a block with the heading "## variables to be set by the user ##". This script assumes that the covariates are Gaussian fit parameters like depth, width, a and b (these begin with "fit_gauss") or the projected HG coefficients (these begin with "proj_hg_coeff_"). These are only use for naming the model, whereas any custom covariates can be used but will result in the model name being "Gauss=none_HG=none". Feel free to edit this naming convention using lines 55-64 under "# name the current model, using Gaussian fit parameters and hg covariateNames".

> To run the model, go into the terminal at the working directory and run: Rscript rv_lm.R "CSV FILE" "VAR1,VAR2,..." If everything runs smoothly a folder called "models" will be created in your working directory with a subfolder of the model's name and metrics that have to do with the fit model itself as a file called "model.rds".  If you want to run the model without any covariates, leave the variables blank.
	
file: rv_lm_leverages.R <br>


*Minimal working example*

This minimal working example will reproduce the results from the Full Model of the paper.

1. Make sure that you have the following packages downloaded: Matrix, tidyverse, stringr, rhdf5.
2. Download the following files: completeLines.csv, rv_lm.R, readSpectra.R.
3. Put all files into the same working directory. In the file "rv_lm.R" replace the WD_DATA object with this working directory and save the file.
4. Go to terminal at this directory, and run: Rscript rv_lm.R completeLines.csv fit_gauss_a,fit_gauss_b,fit_gauss_depth,fit_gauss_sigmasq,proj_hg_coeff_0,proj_hg_coeff_2,proj_hg_coeff_3,proj_hg_coeff_4,proj_hg_coeff_5,proj_hg_coeff_6,proj_hg_coeff_7,proj_hg_coeff_8,proj_hg_coeff_9,proj_hg_coeff_10
5. If everything runs smoothly a folder called "models" will be created in your working directory with a subfolder of the model's name and metrics that have to do with the fit model itself as a file called "model.rds"!

packages: tidyverse, rhdf5, Matrix, patchwork, collapse, parallel, pbmcapply