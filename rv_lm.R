source("readSpectra.R")
library(Matrix)
library(parallel)
library(pbmcapply)

# capture command line arguments
args = commandArgs(trailingOnly = TRUE)

# get csv file name
csvFileName = args[1]
# get covariateNames from command line 
if (is.na(args[2])) {
  covariateNames = c()
} else {
  covariateNames = unlist(strsplit(args[2], ","))
}

## variables to be set by the user ##
# working directory with the data
WD_DATA = "/Users/josephsalzer/research/exostat/"
# response variable
RESPONSE = "rv_template_0.5"
# names of the timeID, lineID, and timeGroupID
TIME_ID_NAME = "date"
LINE_ID_NAME = "line_order"
TIMEGROUP_ID_NAME = "date_groups"
# csv file name
completeLines_df = read_csv(str_c(WD_DATA,csvFileName)) %>%
  mutate(!!TIME_ID_NAME := as.Date( !!sym(TIME_ID_NAME) ))
## variables to be set by the user ##

## get list of lineIDs and timeIDs

# vec of lineIDs  in completeLines_df
lineIDs = completeLines_df %>% group_by(!!sym(LINE_ID_NAME)) %>% summarize(n. = n()) %>% pull(!!sym(LINE_ID_NAME))
# vec of timeIDs in completeLines_df
timeIDs = completeLines_df %>% group_by(!!sym(TIME_ID_NAME)) %>% summarize(n. = n()) %>% pull(!!sym(TIME_ID_NAME))
T_ = length(timeIDs)
L_ = length(lineIDs)

## create model formula and subdirectories

# create directory called "models" in working directory
if (!dir.exists(str_c(WD_DATA,"models"))) {dir.create(str_c(WD_DATA,"models"))}
# base TWFE formula
twfe_formula = paste0("~ ",TIME_ID_NAME," + ",LINE_ID_NAME)
# if covariateNames aren't empty, include interactions, else use TWFE
if (!is_empty(covariateNames)) {
  # build interactions with LINE_ID_NAME for each provided centered covariate
  covar_slopes = paste0(LINE_ID_NAME,":", covariateNames, "_centered", collapse = " + ")
  # construct the full formula
  modelFormula = as.formula(paste(twfe_formula, covar_slopes, sep = " + "))
  
  # name the current model, using Gaussian fit parameters and hg covariateNames
  ## remove or edit these lines if using other covariateNames ##
  gaussCovars = covariateNames[startsWith(covariateNames,"fit_gauss")]
  gaussCovars_names = paste(str_split_i(gaussCovars, "_", i = 3),collapse=",")
  if (gaussCovars_names=="") {gaussCovars_names = "none"} else if(length(gaussCovars) == 4) {gaussCovars_names = "all"}
  hgCovars = covariateNames[startsWith(covariateNames,"proj_hg_coeff_")]
  hgCovars_names = paste(str_split_i(hgCovars, "_", i = 4),collapse=",")
  if (hgCovars_names=="") {hgCovars_names = "none"} else if(length(hgCovars) == 10) {hgCovars_names = "all"}
  model_name = str_c("Gauss=",gaussCovars_names,"_HG=",hgCovars_names)
  
  # create model directory
  model_dir = str_c(WD_DATA,"models/",model_name)
  if (!dir.exists(model_dir)) {dir.create(model_dir)}
  
  rm(twfe_formula,covar_slopes,gaussCovars,gaussCovars_names,hgCovars,hgCovars_names)
} else {
  # if empty, just fit TWFE
  modelFormula = as.formula(twfe_formula)
  # model name
  model_name = "TWFE"
  
  # create model directory
  model_dir = str_c(WD_DATA,"models/",model_name)
  if (!dir.exists(model_dir)) {dir.create(model_dir)}
  rm(twfe_formula)
}

cat("The model formula is:\n", as.character(modelFormula), "\n")
cat("The model name is:", model_name, "\n")

## standardize dataframe

# get standardized rv arranging by line_order and then date
rv_df = ( completeLines_df %>%
            standardize(covariates = covariateNames, response = RESPONSE, lineIDname = LINE_ID_NAME) )$train.df %>%
  arrange(!!sym(LINE_ID_NAME), as.Date(!!sym(TIME_ID_NAME)))

## set contrasts and make design matrix, custom contrasts for timeIDs, sum2zero for lineIDs (single intercept, single slope model)

# create factors out of date + line_order
temp_df = rv_df %>%
  mutate(!!TIME_ID_NAME := factor( !!sym(TIME_ID_NAME) ),
         !!LINE_ID_NAME := factor( !!sym(LINE_ID_NAME) ))
# get group sizes
timeGroup_ids = temp_df %>%
  dplyr::select(!!sym(TIMEGROUP_ID_NAME), !!sym(TIME_ID_NAME)) %>%
  unique()
group_sizes = table(timeGroup_ids[[TIMEGROUP_ID_NAME]])
# (sparse) design matrix for model, set contrasts
designMat = sparse.model.matrix(
  modelFormula,
  temp_df,
  contrasts = setNames(
    list(contr_groupSum(group_sizes), contr.sum),
    c(TIME_ID_NAME, LINE_ID_NAME)
  )
)  
# get every date in each group used to fit the model, this excludes the last date in each group (which should be sum2zero encoded)
timeFE_columns = levels(temp_df[[TIME_ID_NAME]])[-cumsum(group_sizes)]
# append the group ID to the end of these columns
timeFE_columns = str_c(timeFE_columns, "_group", rep( 1:length(group_sizes), group_sizes-1 ) )
# rename columns of design matrix, assuming that the first T are the intercept plus timeFE
# colnames(designMat)[1:sum(group_sizes)] = c("(Intercept)",
#                                             str_c(TIMEGROUP_ID_NAME,seq(1, length(group_sizes)-1)),
#                                             timeFE_columns)
colnames(designMat)[1:sum(group_sizes)] = c("(Intercept)",
                                            if (length(group_sizes) > 1) str_c(TIMEGROUP_ID_NAME,seq(1, length(group_sizes)-1)) else NULL,
                                            timeFE_columns)

# responses
responses = temp_df[[RESPONSE]]

## fit model

# fit the linear model using all data
fit_lm = sparseLM(designMat, responses)

## get time FE results

# initialize 0's in the linear operator
linear_op_mat = Matrix(0, nrow = T_, ncol = length(fit_lm$beta_hat[,1]), sparse = T )
# matrix for estimating the cleaned RV
linear_op_mat[,(length(group_sizes)+1):sum(group_sizes)] = contr_groupSum(group_sizes)[,length(group_sizes):(sum(group_sizes)-1)]
# covariance matrix of model parameters
cov_mat = fit_lm$var_beta_hat
# dataframe of date's cleaned rv
cleanRV_df = data.frame(
  timeID = as.Date(timeIDs)
)
# find the mle estimate, var, and se for the intercept plus alpha
cleanRV_df$estimate = (linear_op_mat %*% fit_lm$beta_hat[,1] )[,1]
cleanRV_covar = linear_op_mat %*% cov_mat %*% t(linear_op_mat)
cleanRV_df$var = diag( cleanRV_covar )
cleanRV_df$se = sqrt(cleanRV_df$var)
# find rmse of the clean RV, comparing the cleaned RV to zero
RMSE = rmse_t(0,cleanRV_df$estimate)

saveRDS(list(designMat = designMat,
             responses = responses,
             df = temp_df,
             timeGroup_ids = timeGroup_ids,
             modelFormula = modelFormula,
             covariateNames = covariateNames,
             fit_lm = fit_lm,
             cleanRV_df = cleanRV_df,
             RMSE = RMSE),
        file = str_c(model_dir, "/model.rds" ) )