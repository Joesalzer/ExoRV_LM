file: readSpectra.R
	r file that contains helper functions for use throughout this project

file: clean_data.Rmd
	Uses the directory line_property_files to join all of the absorption lines together. The output is completeLines.csv, which is a data file for every line's shape properties, line-specific information, and radial velocity. It also provides some summary information and checks of individual lines (like lines that have constant b1 and width, and the lines with smallest rv).
	
file: eda.Rmd
	EDA for this project, includes visuals before models were fit
	