file: readSpectra.R
	r script that contains helper functions for use throughout this project

file: clean_data.Rmd
	Uses the directory "line_property_files" to join all of the absorption lines together. The output is completeLines.csv, which is a data file for every line's shape properties, line-specific information, and radial velocity. It also provides some summary information and checks of individual lines (like lines that have constant b1 and width).
	
file: eda.Rmd
	EDA for this project, includes visuals before models were fit.

file: spectrum_loading.Rmd
	This script includes the opening of various of the raw_spectrum_files. Further, there are tests of the projected hg coefficients and bisectors.

file: modeling.Rmd
	testing grounds for all of our modeling, including the full-data analysis, bootstrap, and cv

file: results.Rmd
	Provide some key results from the full-data analysis, bootstrap, and cv along with visuals for after the modeling was done
		
file: fullData_lm.R
	Fit a two-way fixed effects linear model, with options to include covariates and which csv file to be used. To run the model, go into the terminal at the working directory and run: Rscript rv_lm.R "CSV FILE" "VAR1,VAR2,..." If everything runs smoothly a folder called "models" will be created in your working directory with a subfolder of the model's name and a bunch of things that have to do with the fit model itself as a file called "model.rds".  If you want to run the TWFE model leave the variables blank.