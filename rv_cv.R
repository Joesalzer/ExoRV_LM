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
print( str_c("Number of days in validation set (not including central day): ", 2*CV_NUM_DAYS) )

## USER DEFINED VARS ##
# working directory with the data
WD_DATA = "/Users/josephsalzer/research/exostat/"
# RESPONSE variable
RESPONSE = "rv_template_0.5"
# names of the timeID, lineID, and timeGroupID
TIME_ID_NAME = "date"
LINE_ID_NAME = "line_order"
TIMEGROUP_ID_NAME = "date_groups"
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
if ( !dir.exists(str_c(WD_DATA,"models/",MODEL_NAME,"/cv_block=",2*CV_NUM_DAYS)) ) {
  dir.create(str_c(WD_DATA,"models/",MODEL_NAME,"/cv_block=",2*CV_NUM_DAYS))
}

# list of timeIDs to be left out
left_out_timeIDs = create_sliding_windows(timeIDs, CV_NUM_DAYS)

## IRLS algorithm ##
IRLS = function(rse_tol = 2e-05,
                max.iter = 10,
                X = designMat,
                Y = responses,
                df = completeLines_df,
                lineIDname = LINE_ID_NAME) {
  ## initialize equal weights ##
  iter = 1
  # weight vector
  w_irls = rep(1, nrow(X))
  wls_fit = sparseWLM(X, Y, w_irls)
  # rse for each iteration
  rse_iter = list()
  rse_iter[[1]] = 0
  rse_iter[[2]] = wls_fit$RSE
  
  while( (iter < max.iter) & ( abs(rse_iter[[iter+1]] - rse_iter[[iter]])>rse_tol) ) {
    iter = iter + 1
    # re-weight line-by-line, update weight vector
    w_irls = df %>%
      mutate(r = wls_fit$resid[,1]) %>%
      group_by(!!sym(lineIDname)) %>%
      mutate(sar = sum(abs(r)),
             w = 1/sar) %>%
      ungroup() %>%
      pull(w)
    w_irls = w_irls/sum(w_irls)
    # fit new model
    wls_fit = sparseWLM(X, Y, w_irls)
    rse_iter[[iter + 1]] = wls_fit$RSE
    #print(iter)
  }
  cat("iter: ", iter, " diff RSE: ", rse_iter[[iter+1]] - rse_iter[[iter]], "\n")
  return(list(fit_lm = wls_fit,
              weights = w_irls) )
}

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
  group_sizes_train = train_df %>%
    group_by(!!sym(TIMEGROUP_ID_NAME)) %>%
    summarize(size = n()/L_) %>%
    pull(size)
  group_sizes = df %>%
    group_by(!!sym(TIMEGROUP_ID_NAME)) %>%
    summarize(size = n()/L_) %>%
    pull(size)
  
  # time fixed effects encoding matrix
  if (length(group_sizes_train) > 1) {
    S_time = cbind(
      bdiag(lapply(group_sizes_train, function(n) rep(1, n))) %*% contr.sum(length(group_sizes_train)),
      bdiag(lapply(group_sizes_train, function(n) contr.sum(n) ))
    )
  } else {
    S_time = contr.sum(group_sizes_train)
  }
  # time fixed effects design matrix
  X_train_time = kronecker(rep(1,L_), S_time)
  # design matrix for the line fixed effects
  X_train_line = kronecker( contr.sum(L_, sparse = T), rep(1,T_-length(test_set)) )
  # design matrix for the covariates
  X_train_covar = sparse.model.matrix(modelFormula,train_df) 
  # create full design matrix
  X_train = cbind(rep(1,nrow(X_train_time)),X_train_time,X_train_line,X_train_covar)
  
  if ( startsWith(MODEL_NAME,"LASSO")  ) {
    # get names of columns that are remvoed by lasso
    lasso_removed_names = names(model_fit$lasso_zero_columns[model_fit$lasso_zero_columns])
    # remove LASSO columns
    X_train = X_train[, !(colnames(X_train) %in% lasso_removed_names)]
  }
  rm(X_train_time, X_train_line, X_train_covar)
  
  # responses
  Y_train = train_df[[RESPONSE]]
  
  # fit the linear model using all data, set an appropriate tolerance and maximum number of iterations here
  irls_results = IRLS(rse_tol = 1e-04, max.iter = 10,
                      X = X_train, Y = Y_train, df = train_df)
  fit_lm = irls_results$fit_lm
  train_weights = train_df %>%
    mutate(w = irls_results$weights) %>%
    group_by(!!sym(LINE_ID_NAME)) %>%
    summarize(weight = mean(w))
  #train_weights = irls_results$weights
  
  # get the group-offset encoding for the test set test set based on time point
  g_t = (bdiag(lapply(group_sizes, function(n) rep(1, n))) %*% contr.sum(length(group_sizes)))[(timeIDs %in% test_set),]
  S_test_time = cbind(g_t, Matrix(0, nrow = length(test_set), ncol=sum(group_sizes_train) - length(group_sizes_train) ) )
  X_test_time = kronecker(rep(1,L_), S_test_time)
  
  # design matrix for the line fixed effects
  X_test_line = kronecker( contr.sum(L_, sparse = T), rep(1,length(test_set)) )
  # design matrix for the covariates
  X_test_covar = sparse.model.matrix(modelFormula,test_df) 
  # create full design matrix
  X_test = cbind(rep(1,nrow(X_test_covar)),X_test_time,X_test_line,X_test_covar)
  # responses
  Y_test = test_df[[RESPONSE]]
  
  if ( startsWith(MODEL_NAME,"LASSO")  ) {
    # remove LASSO columns
    X_test = X_test[, !(colnames(X_test) %in% lasso_removed_names)]
  }
  
  rm(X_test_time, X_test_line, X_test_covar)
  
  # create a column for predicted rvs
  test_df[["pred_rv"]] = (X_test %*% fit_lm$beta_hat)[,1]
  
  
  # store data
  saveRDS(
    list(testDF = test_df %>%
           rename(contam_rv = !!sym(RESPONSE) ) %>%
           select(!!sym(TIME_ID_NAME), !!sym(LINE_ID_NAME), contam_rv, pred_rv),
         model_coefs = fit_lm$beta_hat[,1],
         train_weights = train_weights),
    file = str_c(WD_DATA, "models/", MODEL_NAME, "/cv_block=",2*CV_NUM_DAYS,"/cv_df_", test_day, ".rds" ) 
  )
  
}

#cv_parallel(left_out_timeIDs[[1]])
cv_list <- pbmclapply(left_out_timeIDs, FUN = cv_parallel, mc.cores = 5)