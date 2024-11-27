source("readSpectra.R")
library(Matrix)
library(parallel)
library(pbmcapply)

# capture command line arguments
args = commandArgs(trailingOnly = TRUE)

# model name
MODEL_NAME = args[1]
# number of timeIDs in each leave-one-out cv
CV_NUM_DAYS = as.numeric(args[2])
print(CV_NUM_DAYS)

## USER DEFINED VARS ##
# working directory with the data
WD_DATA = "/Users/josephsalzer/research/exostat/"
# RESPONSE variable
RESPONSE = "rv_template_0.5"
# names of the timeID, lineID, and timeGroupID
TIME_ID_NAME = "date"
LINE_ID_NAME = "line_order"
TIMEGROUP_ID_NAME = "date_groups"
## USER DEFINED VARS ##

# read in the model fit
model_fit = readRDS(str_c(WD_DATA, "models/", MODEL_NAME, "/model.rds" ))
df = model_fit$df
covariateNames = model_fit$covariateNames
modelFormula = model_fit$modelFormula

# vec of LineIDs  in completeLines_df
lineIDs = df %>% group_by(!!sym(LINE_ID_NAME)) %>% summarize(n. = n()) %>% pull(!!sym(LINE_ID_NAME))
# vec of timeIDs in completeLines_df
timeIDs = as.Date(df %>% group_by(!!sym(TIME_ID_NAME)) %>% summarize(n. = n()) %>% pull(!!sym(TIME_ID_NAME)))
T_ = length(timeIDs)
L_ = length(lineIDs)

# create directory called "cv" in working directory
if ( !dir.exists(str_c(WD_DATA,"models/",MODEL_NAME,"/cv")) ) {
  dir.create(str_c(WD_DATA,"models/",MODEL_NAME,"/cv"))
}

# list of timeIDs to be left out
left_out_timeIDs = create_sliding_windows(timeIDs, CV_NUM_DAYS)

## CV in parallel ##
cv_parallel = function(i) {
  # get the test set and which day we evaluate on
  test_set = i$test_set
  test_day = i$test_day
  
  # test index if we leave out a given day/several timeIDs
  cv_index = ( df[[TIME_ID_NAME]] %in% test_set )  
  
  # get standardized rv dataframe, seperate into train and test datasets
  standardized_list = standardize(df, covariates = covariateNames, response = RESPONSE, lineIDname = LINE_ID_NAME, test_index = cv_index)
  
  # get training and testing sets
  train_df = standardized_list$train.df %>%
    mutate(!!TIME_ID_NAME := factor( !!sym(TIME_ID_NAME) ),
           !!LINE_ID_NAME := factor( !!sym(LINE_ID_NAME) )) %>%
    arrange(!!sym(LINE_ID_NAME), as.Date(!!sym(TIME_ID_NAME)))
  
  test_df = standardized_list$test.df %>%
    mutate(!!TIME_ID_NAME := factor(!!sym(TIME_ID_NAME)),
           !!LINE_ID_NAME := factor(!!sym(LINE_ID_NAME))) %>%
    arrange(!!sym(LINE_ID_NAME), as.Date(!!sym(TIME_ID_NAME)))
  
  # get group sizes
  timeGroup_ids = train_df %>%
    dplyr::select(!!sym(TIMEGROUP_ID_NAME), !!sym(TIME_ID_NAME)) %>%
    unique()
  group_sizes = table(timeGroup_ids[[TIMEGROUP_ID_NAME]])
  
  # get which date groups are in the testing set
  test_df = test_df %>%
    mutate(!!TIMEGROUP_ID_NAME := factor(!!sym(TIMEGROUP_ID_NAME), levels = names(group_sizes) ) )
  
  # (sparse) design matrix for model, set contrasts
  X_train = sparse.model.matrix(modelFormula,
                                train_df,
                                contrasts = setNames(
                                  list(contr_groupSum(group_sizes), contr.sum),
                                  c(TIME_ID_NAME, LINE_ID_NAME)
                                )
  )  
  # responses
  Y_train = train_df[[RESPONSE]]
  
  # get every date in each group used to fit the model, this excludes the last date in each group (which should be sum2zero encoded)
  timeFE_columns = levels(train_df[[TIME_ID_NAME]])[-cumsum(group_sizes)]
  # append the group ID to the end of these columns
  timeFE_columns = str_c(timeFE_columns, "_group", rep( 1:length(group_sizes), group_sizes-1 ) )
  # rename columns of design matrix
  colnames(X_train)[1:sum(group_sizes)] = c("(Intercept)",
                                            if (length(group_sizes) > 1) {
                                              str_c(TIMEGROUP_ID_NAME,seq(1, length(group_sizes)-1))
                                            } else { NULL },
                                            timeFE_columns)
  
  print("finished making training set")
  
  # fit the linear model on the train data
  fit_train_lm = sparseLM(X_train, Y_train, PRINT_TIME = F)
  
  if ( is.null(fit_train_lm) ) {
    cat("cv failed, singular XtX\n")
    return()
  }
  
  print("finished fitting model")
  
  # (sparse) design matrix for model, without date
  X_test = sparse.model.matrix(
    update(modelFormula, as.formula(paste("~", TIMEGROUP_ID_NAME, "+ . - ",TIME_ID_NAME))),
    test_df,
    contrasts = setNames(list(contr.sum, contr.sum), c(TIMEGROUP_ID_NAME, LINE_ID_NAME)))  
  # responses
  Y_test = test_df[[RESPONSE]]
  
  # get rid of model coefficients associated with date terms (alpha)
  model_coefs_nodate = fit_train_lm$beta_hat[  !grepl("^\\d{4}-\\d{2}-\\d{2}", rownames(fit_train_lm$beta_hat)), ]
  
  # make sure our design matrix and coefficients are the same
  if (!all( names(model_coefs_nodate) == colnames( X_test ) )) {
    cat("error in conforming matrices\n")
    return()
  }
  
  # create a column for predicted rvs
  test_df[["pred_rv"]] = (X_test %*% model_coefs_nodate)[,1] 
  
  # create a column for predicted rvs, subtracting mean RV intercept term (mu)
  #test_df[["pred_rv"]] = (X_test %*% model_coefs_nodate)[,1] - model_coefs_nodate[1]
  
  # store data
  saveRDS(
    list(testDF = test_df %>%
           rename(contam_rv = !!sym(RESPONSE) ) %>%
           select(!!sym(TIME_ID_NAME), !!sym(LINE_ID_NAME), contam_rv, pred_rv),
         model_coefs = model_coefs_nodate),
    file = str_c(WD_DATA, "models/", MODEL_NAME, "/cv/cv_df_", test_day, ".rds" ) 
  )
  
}

cv_list <- pbmclapply(left_out_timeIDs, FUN = cv_parallel, mc.cores = 7)