file: readSpectra.R
	r script that contains helper functions for use throughout this project

file: clean_data.Rmd
	Uses the directory "line_property_files" to join all of the absorption lines' data together. The output is completeLines.csv, which is a data file for every line's shape properties, line-specific information, and radial velocity. It also provides some summary information and checks of individual lines (like lines that have constant b1 and width). Some lines and days are removed from the analysis (and some lines failed to be read-in) as specified in the paper. Projected HG coefficients are also calculated, reading in the file project_template/project_template_deriv_onto_gh.h5
	
file: eda.Rmd
	EDA for this project, includes visuals before models were fit. This file creates some of the figures from the paper: the heatmaps (figure 1) and the RV/depth of select lines (figure 2)

file: spectrum_loading.Rmd
	This script includes the opening of various of the raw_spectrum_files. Further, there are tests of the projected hg coefficients and bisectors.
	
file: pcaRV.Rmd
	This file does PCA on the RV time series of lines.

file: modeling.Rmd
	testing grounds for all of our modeling, including the full-data analysis, bootstrap, and cv

file: results.Rmd
	Provides some key results from the full-data analysis, bootstrap, and cv along with visuals for after the modeling was done
		
file: rv_lm.R
	Fit a two-way fixed effects linear model, with options to include covariates and which csv file to be used. 
	This script has some variables that need to be set by the user (such as the working directory, what the response column is called, what the time column is called, etc.) inside a block with the heading "## variables to be set by the user ##". This script assumes that the covariates are Gaussian fit parameters like depth, width, a and b (these begin with "fit_gauss") or the projected HG coefficients (these begin with "proj_hg_coeff_"). These are only use for naming the model, whereas any custom covariates can be used but will result in the model name being "Gauss=none_HG=none". Feel free to edit this naming convention using lines 54-62 under "# name the current model, using Gaussian fit parameters and hg covariateNames".
	To run the model, go into the terminal at the working directory and run: Rscript rv_lm.R "CSV FILE" "VAR1,VAR2,..." If everything runs smoothly a folder called "models" will be created in your working directory with a subfolder of the model's name and a bunch of things that have to do with the fit model itself as a file called "model.rds".  If you want to run the model without any covariates, leave the variables blank.