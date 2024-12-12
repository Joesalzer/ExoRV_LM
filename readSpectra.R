###########################
# helper functions, to read in NEID datafiles and fit/evaluate models
###########################

require(tidyverse)
require(stringr)
require(rhdf5)
require(Matrix)

###########################
# center, scale, and change function
###########################

centerFun = function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm))
scaleFun = function(x, na.rm = FALSE) (x/sd(x, na.rm = na.rm))
changeFun = function(x) ( x - x[1] )

###########################
# load a "solar_spectrum_DATE.h5" file into a dataframe #
###########################

# spectra_name: string for file name
# wd_data: working directory

loadSolarFile = function(spectra_name, wd_data = "/Users/josephsalzer/research/exostat/raw_spectrum_files/") {
  
  # download solar spectrum, with flux, flux_norm, lambda, var, and var_norm
  # (pixels, order_idx)
  solar_spectrum = H5Fopen(str_c(wd_data, spectra_name))
  
  # number of pixels and orders
  num_pixels = dim(solar_spectrum$flux)[1]
  num_orders = dim(solar_spectrum$flux)[2]
  
  # vector of orders and pixels
  orders = seq(1,num_orders)
  pixels = seq(1,num_pixels)
  
  # create dataframe for orders and pixels
  solar_spectrum_df = expand_grid(order = orders, pixel = pixels)
  # set columns
  solar_spectrum_df$flux_norm = as.vector(solar_spectrum$flux_norm)
  solar_spectrum_df$flux = as.vector(solar_spectrum$flux)
  solar_spectrum_df$lambda = as.vector(solar_spectrum$lambda)
  solar_spectrum_df$var = as.vector(solar_spectrum$var)
  solar_spectrum_df$var_norm = as.vector(solar_spectrum$var_norm)
  solar_spectrum_df$date = as.Date( str_split_1(spectra_name,"_")[3] )
  
  return(solar_spectrum_df)
  
}

###########################
## load a "line_LINE_ORDER.h5: file into a dataframe ##
###########################

# line_name: string for file name
# wd_data: working directory

loadLineFile = function(line_name, wd_data = "/Users/josephsalzer/research/exostat/raw_spectrum_files/") {
  # download a single line
  # (pixels,date)
  single_line = H5Fopen(str_c(wd_data, line_name))
  
  # vector of dates
  dates = as.Date(single_line$dates)
  # vector of pixels
  pixels = as.numeric( str_split_1(single_line$pixels,":") )
  pixels = seq(pixels[1],pixels[2])
  
  # create dataframe for dates and pixels
  single_line.df = expand_grid(date = dates, pixel = pixels )
  # set columns
  single_line.df$flux = as.vector( single_line$flux)
  single_line.df$flux_norm = as.vector( single_line$flux_norm)
  single_line.df$lambda = as.vector( single_line$lambda) 
  single_line.df$var = as.vector( single_line$var)
  single_line.df$var_norm = as.vector( single_line$var_norm)
  single_line.df$lambdac_fit = single_line$lambdac_fit
  single_line.df$lambdac_vald = single_line$lambdac_vald
  #single_line.df$pixels = single_line$pixels
  single_line.df$order_idx = single_line$order_idx
  single_line.df$prder_phys = single_line$order_phys
  single_line.df$deltav = single_line$deltav
  
  return(single_line.df)
}

###########################
## load a line_property_file into a dataframe ##
###########################

# file returns a list with
# line.df: contains the line-day-measurements as a dataframe
# template.df: contains the template of the line's flux and wavelength
# reconstruct.df: contains the gh/tgh reconstructed line for each day
# line_name: string for file name
# wd_data: working directory

loadLineShape = function(line_name, wd_data = "/Users/josephsalzer/research/exostat/line_property_files/") {
  
  # open one file
  line.h5 = H5Fopen(str_c(wd_data, line_name))
  # see names, datatype and dim of line.h5
  name.df = h5ls(line.h5)
  # list of values for each name in line.h5
  value.list = h5dump(line.h5)
  
  # names from the name.df to include in line.df, only the 345 x 1 
  namesDate = name.df$name[name.df$dim == "345"]
  # names from the name.df to include in line.df, only the constant values
  namesConstant = name.df$name[name.df$dim == "( 0 )"]
  # names from the name.df to include in line.df, only the size of the number of HG x day (11 x 345) included
  namesHGxDate = name.df$name[name.df$dim ==  "11 x 345"]
  
  # go to next if the there aren't 4 names with HG by date
  if ( length( namesHGxDate ) != 4 ) {
    print("FAILED TO READ")
    return(list("line.df" = NA))
  }
  
  # create line.df file of shape measurements/rv over time
  line.df = data.frame( value.list[namesDate] )
  
  # create hg dataframe, binding the namesHGxDate into a 345 x 44 dataframe
  hg.df = data.frame( t( do.call(rbind, value.list[namesHGxDate]) ) )
  # rename columns
  colnames(hg.df) = paste( rep(gsub("gh", "hg", namesHGxDate), each = 11), "_", c(0:10), sep = "")
  
  # bind hg.df, line.df, and the constants 
  line.df = cbind(line.df, hg.df, do.call(cbind, value.list[namesConstant] ))
  
  # add line-order column to line.df
  line.df$line_order = str_flatten(str_split_1( line_name, "_" )[2:3], collapse = "_")
  
  # replace σ with sigma in colnames
  colnames(line.df) = gsub("σ","sigma_",colnames(line.df))
  
  # get template line dataframe
  template.df = data.frame(
    flux = value.list$template_flux,
    lambda = value.list$template_lambda,
    var = value.list$template_var
  )
  
  # get reconstructed line dataframe
  reconstruct.df = data.frame(
    date = value.list$date,
    gh_reconstruct = t(value.list$gh_reconstruct),
    tgh_reconstruct = t(value.list$tgh_reconstruct)
  )
  
  return( list("line.df" = line.df,
               "template.df" = template.df,
               "reconstruct.df" = reconstruct.df) )
}

###########################
## Standardize our data and train-test split ##
###########################

# input:
# current_df - dataframe with lineID, date, and various covariates and a response as columns
# covariates - vector of strings for covariates that we want to standardize
# test_index - vector of logical for the indices that we would like to split from training data
# response - string for the response column
# lineIDname - string for the lineID

standardize = function(current_df, covariates = NULL, test_index = NULL, response = "rv_template_0.5", lineIDname = "line_order") {
  
  # if we want to do training and testing, include indices in test_index and this makes a train and test set,
  # otherwise it just makes a train set of the current df
  if (!is.null( test_index )) {
    # create training data without test_index
    train = current_df[!test_index,]
    # create testing data with test_index
    test = current_df[test_index,]
  } else {
    train = current_df
  }
  
  # if no covariates are included, return original df
  if (is.null(covariates)) {
    if (is.null( test_index )) {
      return( list(train.df = train) )
    } else {
      return( list(train.df = train,
                   test.df = test) )
    }
  }
  
  ## this section normalizes train data by means and sd of train data ##
  
  # get means of our covariates per line (for training days)
  train_means = train %>%
    group_by(!!sym(lineIDname)) %>%
    summarize_at(covariates, mean)%>%
    rename_at(covariates, ~ paste0(., '_mean'))
  # get sds of our centered covariates (for training days)
  train_sds = train %>%
    group_by(!!sym(lineIDname)) %>%
    mutate_at(covariates,centerFun) %>%
    ungroup() %>%
    summarize_at(covariates,sd)
  # merge the training shape measurement means by lineIDname
  train = merge(train, train_means, by = lineIDname)
  # loop over covariates, centering by the mean
  for (covar in covariates) {
    train[[paste0(covar, "_centered")]] = train[[covar]] - train[[paste0(covar, "_mean")]]
  }
  # divide by standard deviation of train data
  train[,str_c(covariates,"_centered")] = train[,str_c(covariates,"_centered")]/train_sds[rep(1,dim(train)[1]),]
  # drop mean covariates from train data
  train = train[, !colnames(train) %in% str_c(covariates,"_mean") ]

  ## if we dont want to do training and testing, this just returns the train set ##
  if (is.null( test_index )) {
    return(list(train.df = train,
                train_means = train_means,
                train_sds = train_sds))
  }

  ## this section normalizes test data by means and sd of train data ##
  
  # merge the training shape measurement means by lineIDname
  test = merge(test, train_means, by = lineIDname)
  # loop over covariates, centering by the mean
  for (covar in covariates) {
    test[[paste0(covar, "_centered")]] = test[[covar]] - test[[paste0(covar, "_mean")]]
  }
  # divide by standard deviation of train data
  test[,str_c(covariates,"_centered")] = test[,str_c(covariates,"_centered")]/train_sds[rep(1,dim(test)[1]),]
  # drop mean covariates from test data
  test = test[, !colnames(test) %in% str_c(covariates,"_mean") ]
  
  return(list(train.df = train,
              test.df = test,
              train_means = train_means,
              train_sds = train_sds))
}


###########################
## fit a LM using model matrix and responses##
###########################

# X (n x p) : model matrix
# Y (n x 1) : responses
# PRINT_TIME :logical for whether we print of not

# returns a list of model fits
sparseLM = function(X, Y, PRINT_TIME = T) {
  START_TIME = Sys.time()
  
  # get XtX
  XtX = crossprod(X)
  
  # get X^T Y
  XtY = crossprod(X, Y)
  
  # get (X^T X)^-1 matrix
  # ...using cholesky decomp
  XtX_inv =  try( chol2inv(Matrix::Cholesky(XtX)) )
  
  # catch error if XtX is singular
  if (all(class(XtX_inv) == "try-error") ) {
    cat("singular XtX\n")
    return()
  }
  
  # beta hat, linear model coefficients
  beta_hat = XtX_inv %*% XtY
  
  # fitted values and residuals
  y_hat = X %*% beta_hat
  resid = Y - y_hat
  
  # dimensions
  p = ncol(X)
  n = nrow(X)
  
  # sigma2 hat and variance of beta_hat
  sigma2_hat = drop( crossprod(resid) / (n-p) )
  var_beta_hat = sigma2_hat * XtX_inv
  
  # sum of squared residuals
  SSR = sum(resid^2)
  # total sum of squares
  SST = sum( ( Y - mean(Y) )^2 )
  # R-squared, AIC, BIC calculations
  R2 = 1-SSR/SST
  adj_R2 = 1 - (1- R2)*(n - 1)/(n - p - 1)
  AIC = n * log(SSR/n) + 2 * p
  BIC = n * log(SSR/n) + 2 * log(n)
  # residual standard errors
  RSE = sqrt( SSR/(n-p) )
  
  if (PRINT_TIME) {
    print(Sys.time() - START_TIME)
  }
  

  return( list( beta_hat = beta_hat,
                y_hat = y_hat,
                resid = resid,
                sigma2_hat = sigma2_hat,
                var_beta_hat = var_beta_hat,
                adj_R2 = adj_R2,
                AIC = AIC,
                BIC = BIC,
                RSE = RSE))
}

###########################
## case resampling for model coefficients see Algorithm 6.2 in Davison and Hinkley ##
###########################

# runs a single instance of the case re sampling bootstrap

# i : an int that sets the seed for replication
# X (n x p) : model matrix, defaults to "designMat"
# Y (n x 1) : responses, defaults to "responses"
# wd : string of working directory, defaults to wd_data
# model : string of model name, defaults ot model_name

# returns an R x p matrix of the bootstrap coefficients
casebootCoef_parallel_old = function(i, X = designMat, Y = responses, wd = wd_data, model = model_name) {
  set.seed(i)
  # sample with replacement of indices 
  idx_boot = sample( seq(1,nrow(X)) , replace = T)
  
  # boot strap the design matrix and response
  X_boot = X[idx_boot,]
  Y_boot = Y[idx_boot]
  
  # fit the bootstrap matrices
  fit_boot = sparseLM(X_boot, Y_boot, PRINT_TIME = F)
  # store bootstrap coefficients
  coef_boot = fit_boot$beta_hat[,1]
  # store bootstrap residuals
  resid_boot = fit_boot$resid[,1]
  
  # store data
  saveRDS(list(coef_boot = coef_boot),
          file = str_c(wd, "models/", model, "/bootstraps/bootstraps_list_", i, ".rds" ) )
  
  # return( list(coef_boot = coef_boot,
  #              resid_boot = resid_boot,
  #              indices_boot = idx_boot) )
  
}

###########################
## wild bootstrap for model coefficients ##
###########################

# df : already-standardized dataframe  with timeIDname and lineIDname
# fittedLM : already fit linear model on the dataframe
# X : design matrix
# name : name of the model
# timeIDname : column name of the time variable
# lineIDname : column name of the lines
# resid : residuals of the linear model, can be normalized by leverage beforehand
# wd : working directory to store bootstrapped coefficients

# returns an R x p matrix of the bootstrap coefficients
wildBoot_parallel_old = function(seedID, df = rv_df, fittedLM = fit_lm, X = designMat, name = model_name,
                             timeIDname = "date", lineIDname = "line_order", resid = u_hat, wd = "/Users/josephsalzer/research/exostat/") {
  
  # set a given seed for replication
  set.seed(seedID)
  
  # ensure that the df is arranged first by lineID and then timeID, and that these are factors
  df = df %>%
    arrange(!!sym(lineIDname), !!sym(timeIDname)) %>%
    mutate( lineID = factor(!!sym(lineIDname)),
            timeID = factor(!!sym(timeIDname)))
  
  # rademacher random variable for each block in the dataset
  v = sample( c(-1, 1), size = length(unique(df$boot_block)), replace = TRUE)
  # generate wild bootstrap residuals by block
  u_star = resid*v[df$boot_block]
  
  # construct bootstrap sample
  y_star = fittedLM$y_hat + u_star
  # re-estimate the model on bootstrap sample
  fitstar_lm = sparseLM(X, y_star, PRINT_TIME = F)

  # store data
  saveRDS(coef_boot = fitstar_lm$beta_hat[,1],
          file = str_c(wd, "models/", name, "/wild_bootstraps/bootstraps_list_", seedID, ".rds" ) )
  
}

###########################
# function to create list of of days to do cv over, a sliding window over the days (not based on indices but the days themselves), these are overlapping
###########################

# dates : vector of days to do cv on sequence
# window_size : number of days to use in each CV

create_sliding_windows = function(dates, window_size) {
  # initialize list of validations sets
  val_sets = list()
  
  # loop through all days and create windows
  for (day in dates) {
    # get the days in our dataset that are within that window
    val_set = list(test_set = dates[(dates >= day - window_size) & (dates <= day + window_size)],
                   test_day = as.Date(day))
    
    
    val_sets = append(val_sets, list(val_set) )
  }
  
  return(val_sets)
}

###########################
## run a single instance of the leave-one day out cv ##
###########################

# i : the days that we want to leave out of the model
# df : dataframe to be passed to several functions, requires a line_order, date, rv_template_0.5, and other covars
# days : vector of days in the dataframe
# covariates : vector of covariates in the model
# covariates_pca : vector of pca covariates in the model
# model_form : the model formula to use
# wd : string of working directory, defaults to wd_data
# model : string of model name, defaults ot model_name

# returns a dataframe of date, line_order, contam_rv, and pred_rv
cv_parallel_old = function(i, df = completeLines_df,
                       covariates = covariateNames,
                       model_form = modelFormula, wd = wd_data, model = model_name,
                       response = "rv_template_0.5",
                       timeIDname = "date",
                       lineIDname = "line_order",
                       timeGroupIDname = "date_groups") {
  # get the test set and which day we evaluate on
  test_set = i$test_set
  test_day = i$test_day
  
  # test index if we leave out a given day/several timeIDs
  cv_index = ( df[[timeIDname]] %in% test_set )  
  
  # get standardized rv dataframe, seperate into train and test datasets
  standardized_list = standardize(df, covariates = covariateNames, test_index = cv_index)
  
  # get training and testing sets
  train_df = standardized_list$train.df %>%
    mutate(!!timeIDname := factor( !!sym(timeIDname) ),
           !!lineIDname := factor( !!sym(lineIDname) )) %>%
    arrange(!!sym(lineIDname), as.Date(!!sym(timeIDname)))
  
  test_df = standardized_list$test.df %>%
    mutate(!!timeIDname := factor(!!sym(timeIDname)),
           !!lineIDname := factor(!!sym(lineIDname))) %>%
    arrange(!!sym(lineIDname), as.Date(!!sym(timeIDname)))
  
  # get group sizes
  # group_sizes = train_df %>%
  #   dplyr::select(!!sym(timeGroupIDname), !!sym(timeIDname)) %>%
  #   unique() %>%
  #   group_by(!!sym(timeGroupIDname)) %>%
  #   summarise(n = n()) %>%
  #   pull(n)
  # get group sizes
  timeGroup_ids = train_df %>%
    dplyr::select(!!sym(timeGroupIDname), !!sym(timeIDname)) %>%
    unique()
  group_sizes = table(timeGroup_ids[[timeGroupIDname]])
  
  # get which date groups are in the testing set
  # test_df = test_df %>%
  #   mutate(!!timeGroupIDname := factor(!!sym(timeGroupIDname), levels = c(1:length(group_sizes)) ) )
  test_df = test_df %>%
    mutate(!!timeGroupIDname := factor(!!sym(timeGroupIDname), levels = names(group_sizes) ) )
  
  # (sparse) design matrix for model, set contrasts
  X_train = sparse.model.matrix(model_form,
                                train_df,
                                contrasts = setNames(
                                  list(contr_groupSum(group_sizes), contr.sum),
                                  c(timeIDname, lineIDname)
                                )
  )  
  # responses
  Y_train = train_df[[response]]
  
  # get every date in each group used to fit the model, this excludes the last date in each group (which should be sum2zero encoded)
  timeFE_columns = levels(train_df[[timeIDname]])[-cumsum(group_sizes)]
  # append the group ID to the end of these columns
  timeFE_columns = str_c(timeFE_columns, "_group", rep( 1:length(group_sizes), group_sizes-1 ) )
  # rename columns of design matrix
  colnames(X_train)[1:sum(group_sizes)] = c("(Intercept)",str_c(timeGroupIDname,seq(1, length(group_sizes)-1) ),timeFE_columns)
  
  print("finished making training set")
  
  # fit the linear model on the train data
  fit_train_lm = sparseLM(X_train, Y_train, PRINT_TIME = F)
  
  if ( is.null(fit_train_lm) ) {
    cat("cv failed, singular XtX\n")
    return()
  }
  
  print("finished fitting model")
  
  # (sparse) design matrix for model, without date
  X_test = sparse.model.matrix(update(model_form, as.formula(paste("~", timeGroupIDname, "+ . - ",timeIDname))),
                               test_df,
                               contrasts = setNames(list(contr.sum, contr.sum), c(timeGroupIDname, lineIDname)))  
  # responses
  Y_test = test_df[[response]]
  
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
  
  # create directory called "cv" in working directory
  if ( !dir.exists(str_c(wd_data,"models/",model,"/cv")) ) {
    dir.create(str_c(wd_data,"models/",model,"/cv"))
  }
  # store data
  saveRDS(
    list(testDF = test_df %>%
           rename(contam_rv = !!sym(response) ) %>%
           select(!!sym(timeIDname), !!sym(lineIDname), contam_rv, pred_rv),
         model_coefs = model_coefs_nodate),
    file = str_c(wd, "models/", model, "/cv/cv_df_", test_day, ".rds" ) 
  )
  
}

###########################
## inject a planet of amplitude and frequency into our dataframe ##
###########################

# df : datframe that includes line_order, date, and rv_template_0.5
# amp : amplitude of the sin wave (defaults to 1)
# freq : frequency of rotation (defaults to earth-like rotation of 366 days)
# offset : horizontal shift of the planet signal (defaults to none)

injectPlanet = function(df, amp = 1, freq = 2*pi/366, horizontal_offset = 0, vertical_offset = 0, lineIDname = "line_order", timeIDname = "date", response = "rv_template_0.5") {
  df %>%
    # for every line, get the date as a day value (starting at 0)
    group_by(!!sym(lineIDname)) %>%
    mutate(time_num = changeFun(as.numeric(!!sym(timeIDname))) ) %>%
    ungroup() %>%
    # for every day calculate the amount to perturb a given line
    group_by(!!sym(timeIDname)) %>%
    mutate(planet_val = vertical_offset + amp*sin(freq*time_num + horizontal_offset) ) %>%
    ungroup() %>%
    # add perturbation to each line-day
    mutate(planet_rv = !!sym(response) + planet_val) %>%
    # include mean perturbed rv for each day
    group_by(!!sym(timeIDname)) %>%
    mutate(planet_rv_time = mean(planet_rv)) %>%
    ungroup() %>%
    # include mean perturbed rv for each line
    group_by(!!sym(lineIDname)) %>%
    mutate(planet_rv_line = mean(planet_rv)) %>%
    ungroup()
}

###########################
## calculate RMSE for a vector ##
###########################
# c - offset from
# v - vector to get rmse wrt to the offset
rmse_t = function(c, v) { sqrt( mean( (v - c)^2 ) ) }

# ###########################
# ## function to place a list of matrices on the diagonal ##
# ###########################
# # matrices: list of matrices to put on the diagonal
# block_diagonal = function(matrices) {
#   # calculate the total number of rows and columns for the resulting matrix
#   total_rows = sum(sapply(matrices, nrow))
#   total_cols = sum(sapply(matrices, ncol))
#   # initialize zero matrix
#   result = matrix(0, nrow = total_rows, ncol = total_cols)
#   # fill the diagonal with the input matrices
#   current_row = 1
#   current_col = 1
#   for (mat in matrices) {
#     rows = nrow(mat)
#     cols = ncol(mat)
#     # place the matrix in the current diagonal block
#     result[current_row:(current_row + rows - 1), current_col:(current_col + cols - 1)] = mat
#     # update
#     current_row = current_row + rows
#     current_col = current_col + cols
#   }
#   return(result)
# }

###########################
## function to make a custom contrast for time FE ##
###########################
# this ensures that each treatment is with respect to its own group's mean and they are sum to zero encoded within each group
# group_sizes: vector of group sizes such that each treatment within each group is given a mean offset
contr_groupSum = function(group_sizes) {
  # number of groups
  num_groups = length(group_sizes)
  # if only one group, return sum2zero contrasts for time
  if (num_groups == 1) {
    return(contr.sum(group_sizes))
  }
  # return contrast, the first columns encode the grand means for each group, and the subsequent columns encode the sum2zero contrasts within each group
  return(
    cbind(
      bdiag( lapply(group_sizes, FUN = function(x) rep(1,x) ) ) %*% contr.sum(num_groups),
      bdiag(lapply(group_sizes, contr.sum))
    )
  )
}

# function to make a custom contrast from a vector of group sizes and number of lines
# the resulting contrast is a num_groups*num_lines by num_groups*(num_lines-1) matrix that encodes each line order has its own intercept for each date group
lineGroupFE_contrast = function(group_sizes, num_lines) {
  # number of groups
  num_groups = length(group_sizes)
  # number of time points
  num_times = sum(group_sizes)
  # produce a lineOrder ID matrix, encoding which line a given row is
  lineOrder_ID = matrix( rep(1,num_groups), nrow=1) %x% contr.sum(num_lines) %x% matrix( rep(1,num_groups), ncol=1)
  # produce a group_ID matrix, encodigin which group a given row is in
  group_ID = cbind(
    rep(1,num_groups),
    contr.sum(num_groups)
  )
  group_ID = matrix( rep(1,num_lines), ncol=1) %x% group_ID %x% matrix( rep(1,num_lines-1), nrow=1)
  # return contrast
  return( lineOrder_ID*group_ID )
}

# cv_parallel_old = function(i, df = final_df,
#                                  covariates = covars,
#                                  model_form = modelFormula, wd = wd_data, model = model_name,
#                                  response = "rv_template_0.5", timeIDname = "date", lineIDname = "line_order") {
#   # get the test set and which day we evaluate on
#   test_set = i$test_set
#   test_day = i$test_day
#   
#   # test index if we leave out a given day/several days
#   cv_index = ( df[[timeIDname]] %in% test_set )  
#   
#   # get standardized rv dataframe, seperate into train and test datasets
#   standardized_list = standardize(df, covariates = covariates, test_index = cv_index)
#   
#   # get training and testing sets
#   train_df = standardized_list$train.df %>%
#     mutate(date = factor( !!sym(timeIDname) ),
#            line_order = factor( !!sym(lineIDname) )) %>%
#     arrange(!!sym(lineIDname), as.Date(!!sym(timeIDname)))
#   
#   test_df = standardized_list$test.df%>%
#     mutate(date = factor(!!sym(timeIDname)),
#            line_order = factor(!!sym(lineIDname))) %>%
#     arrange(!!sym(lineIDname), as.Date(!!sym(timeIDname)))
# 
#   # get group sizes
#   group_sizes = train_df %>%
#     dplyr::select(date_groups, !!sym(timeIDname)) %>%
#     unique() %>%
#     group_by(date_groups) %>%
#     summarise(n = n()) %>%
#     pull(n)
#   
#   # (sparse) design matrix for model, set contrasts
#   X_train = sparse.model.matrix(model_form,train_df,
#                                 contrasts = list(date = contr_groupSum(group_sizes),
#                                                    line_order = contr.sum))  
#   # responses
#   Y_train = train_df[[response]]
#   
#   # rename columns of design matrix
#   colnames(X_train)[1:sum(group_sizes)] = c("intercept","groupID",str_c("date",seq(1,group_sizes[1]-1),"_group1"),str_c("date",seq(1,group_sizes[2]-1),"_group2"))
#   
#   print("finished making training set")
#   
#   # fit the linear model on the train data
#   fit_train_lm = sparseLM(X_train, Y_train, PRINT_TIME = F)
#   
#   if ( is.null(fit_train_lm) ) {
#     cat("cv failed, singular XtX\n")
#     return()
#   }
#   
#   print("finished fitting model")
# 
#   # (sparse) design matrix for model, without date
#   # X_test = Matrix(
#   #   model.matrix(update(model_form, ~ . - date), test_df, contrasts.arg = list(line_order = "contr.sum")),
#   #   sparse = T)
#   X_test = sparse.model.matrix(update(model_form, ~ date_groups + . - date ),test_df,contrasts = list(line_order = contr.sum))  
#   # responses
#   Y_test = test_df[[response]]
#   
#   # get rid of model coefficients associated with date terms (alpha)
#   model_coefs_nodate = fit_train_lm$beta_hat[ !startsWith( rownames(fit_train_lm$beta_hat), "date"), ]
#   
#   # make sure our design matrix and coefficients are the same
#   if (!all( names(model_coefs_nodate) == colnames( X_test ) )) {
#     cat("error in conforming matrices\n")
#     return()
#   }
#   
#   # create a column for predicted rvs
#   test_df[["pred_rv"]] = (X_test %*% model_coefs_nodate)[,1] 
#   
#   # create a column for predicted rvs, subtracting mean RV intercept term (mu)
#   #test_df[["pred_rv"]] = (X_test %*% model_coefs_nodate)[,1] - model_coefs_nodate[1]
#   
#   # store data
#   saveRDS(
#     list(testDF = test_df %>%
#            rename(contam_rv = !!sym(response) ) %>%
#            select(!!sym(timeIDname), !!sym(lineIDname), contam_rv, pred_rv),
#          model_coefs = model_coefs_nodate),
#     file = str_c(wd, "models/", model, "/cv/cv_df_", test_day, ".rds" ) 
#   )
#   
# }
# # X (n x p) : model matrix
# # Y (n x 1) : responses
# # w (L x 1) : weights for lines
# # PRINT_TIME :logical for whether we print of not
# 
# # returns a list of model fits
# 
# sparseWLM = function(X, Y, w, PRINT_TIME = T) {
#   START_TIME = Sys.time()
#   # number of dimensions
#   p = dim(X)[2]
#   # number of dimensions
#   n = dim(X)[1]
#   # number of lines
#   L_ = length(w)
#   # number of time points
#   T_ = n/L_
#   
#   # weight matrix
#   W = rep( w, each = T_ )*Diagonal( T_*L_ )
#   
#   # get Xt W X
#   XtWX = t(X) %*% W %*% X
#   
#   # get (X^T W X)^-1 matrix
#   # ...using cholesky decomp
#   XtWX_inv =  try( chol2inv(chol(XtWX)) )
#   
#   # catch error if XtX is singular
#   if (all(class(XtWX_inv) == "try-error") ) {
#     cat("singular XtX\n")
#     return()
#   }
#   
#   # get X^T W Y
#   XtWY = t(X) %*% W %*% Y
#   
#   # beta hat, linear model coefficients
#   beta_hat = XtWX_inv %*% XtWY
#   
#   # fitted value
#   y_hat = X %*% beta_hat
#   
#   # residuals
#   resid = Y - y_hat
#   
#   # sigma2 hat, 
#   sigma2_hat = t(resid) %*% W %*% resid/(length(Y) - p)
#   
#   # covariance matrix of coefficients
#   var_beta_hat = drop(sigma2_hat)*XtWX_inv
#   
#   if (PRINT_TIME) {
#     print( Sys.time() - START_TIME )
#   }
#   
#   # sum of squared residuals
#   SSR = sum(resid^2)
#   # total sum of squares
#   SST = sum( ( Y - mean(Y) )^2 )
#   # r squared
#   R2 = 1-SSR/SST
#   # adjusted r squared
#   adj_R2 = 1 - (1- R2)*(n - 1)/(n - p - 1)
#   # AIC
#   AIC = n * log(SSR/n) + 2 * p
#   # BIC
#   BIC = n * log(SSR/n) + 2 * log(n)
#   
#   
#   return( list( beta_hat = beta_hat,
#                 y_hat = y_hat,
#                 resid = resid,
#                 sigma2_hat = sigma2_hat,
#                 var_beta_hat = var_beta_hat,
#                 adj_R2 = adj_R2,
#                 AIC = AIC,
#                 BIC = BIC,
#                 W = W))
# }
# # returns a dataframe of date, line_order, contam_rv, and pred_rv
# cv_parallel_weighted = function(i, df = final_df,
#                        covariates = covars,
#                        model_form = modelFormula, wd = wd_data, model = model_name,
#                        response = "rv_template_0.5", timeIDname = "date", lineIDname = "line_order",
#                        w = weightByLine) {
#   # get the test set and which day we evaluate on
#   test_set = i$test_set
#   test_day = i$test_day
#   
#   # test index if we leave out a given day/several days
#   cv_index = ( df[[timeIDname]] %in% test_set )  
#   
#   # get standardized rv dataframe, seperate into train and test datasets
#   standardized_list = standardize(df, covariates = covariates, test_index = cv_index)
#   
#   # get training and testing sets
#   train_df = standardized_list$train.df %>%
#     mutate(date = factor( !!sym(timeIDname) ),
#            line_order = factor( !!sym(lineIDname) ) ) %>%
#     arrange(!!sym(lineIDname), !!sym(timeIDname))
#   test_df = standardized_list$test.df%>%
#     mutate(date = factor(!!sym(timeIDname)),
#            line_order = factor(!!sym(lineIDname))) %>%
#     arrange(!!sym(lineIDname), !!sym(timeIDname))
#   
#   # (sparse) design matrix for model, using sum to zero contrasts for days + lines
#   X_train = Matrix(
#     model.matrix(model_form, train_df, contrasts.arg = list(date = "contr.sum", line_order = "contr.sum")),
#     sparse = T)
#   # responses
#   Y_train = train_df[[response]]
#   
#   print("finished making training set")
#   
#   # fit the linear model on the train data
#   fit_train_lm = sparseWLM(X_train, Y_train, w, PRINT_TIME = F)
#   
#   if ( is.null(fit_train_lm) ) {
#     cat("cv failed, singular XtX\n")
#     return()
#   }
#   
#   print("finished fitting model")
#   
#   # (sparse) design matrix for model, without date
#   X_test = Matrix(
#     model.matrix(update(model_form, ~ . - date), test_df, contrasts.arg = list(line_order = "contr.sum")),
#     sparse = T)
#   # responses
#   Y_test = test_df[[response]]
#   
#   # get rid of model coefficients associated with date terms (alpha)
#   model_coefs_nodate = fit_train_lm$beta_hat[ !startsWith( rownames(fit_train_lm$beta_hat), timeIDname), ]
#   
#   # make sure our design matrix and coefficients are the same
#   if (!all( names(model_coefs_nodate) == colnames( X_test ) )) {
#     cat("error in conforming matrices\n")
#     return()
#   }
#   
#   # create a column for predicted rvs
#   test_df[["pred_rv"]] = (X_test %*% model_coefs_nodate)[,1] 
#   
#   # create a column for predicted rvs, subtracting mean RV intercept term (mu)
#   #test_df[["pred_rv"]] = (X_test %*% model_coefs_nodate)[,1] - model_coefs_nodate[1]
#   
#   # store data
#   saveRDS(
#     list(testDF = test_df %>%
#            rename(contam_rv = !!sym(response) ) %>%
#            select(!!sym(timeIDname), !!sym(lineIDname), contam_rv, pred_rv),
#          model_coefs = model_coefs_nodate),
#     file = str_c(wd, "models/", model, "/cv_weighted/cv_df_", test_day, ".rds" ) 
#   )
#   
#   
# }
# sparseLM = function(X, Y, PRINT_TIME = T) {
#   START_TIME = Sys.time()
#   
#   # get XtX
#   XtX = t(X) %*% X
#   
#   # get (X^T X)^-1 matrix
#   # ...using cholesky decomp
#   XtX_inv =  try( chol2inv(chol(XtX)) )
#   
#   # catch error if XtX is singular
#   if (all(class(XtX_inv) == "try-error") ) {
#     cat("singular XtX\n")
#     return()
#   }
#   
#   # get X^T Y
#   XtY = t(X) %*% Y
#   
#   # beta hat, linear model coefficients
#   beta_hat = XtX_inv %*% XtY
#   
#   # fitted value
#   y_hat = X %*% beta_hat
#   
#   # residuals
#   resid = Y - y_hat
#   
#   # number of dimensions
#   p = dim(X)[2]
#   # number of dimensions
#   n = dim(X)[1]
#   
#   # sigma2 hat, 
#   sigma2_hat = t(resid) %*% resid/(length(Y) - p)
#   
#   # covariance matrix of coefficients
#   var_beta_hat = drop(sigma2_hat)*XtX_inv
#   
#   if (PRINT_TIME) {
#     print( Sys.time() - START_TIME )
#   }
#   
#   # sum of squared residuals
#   SSR = sum(resid^2)
#   # total sum of squares
#   SST = sum( ( Y - mean(Y) )^2 )
#   # r squared
#   R2 = 1-SSR/SST
#   # adjusted r squared
#   adj_R2 = 1 - (1- R2)*(n - 1)/(n - p - 1)
#   # AIC
#   AIC = n * log(SSR/n) + 2 * p
#   # BIC
#   BIC = n * log(SSR/n) + 2 * log(n)
#   
#   
#   return( list( beta_hat = beta_hat,
#                 y_hat = y_hat,
#                 resid = resid,
#                 sigma2_hat = sigma2_hat,
#                 var_beta_hat = var_beta_hat,
#                 adj_R2 = adj_R2,
#                 AIC = AIC,
#                 BIC = BIC))
# }
# injectPlanet2 = function(df, amp = 1, freq = 2*pi/366, offset = 0) {
#   df %>%
#     # for every line, get the date as a day value (starting at 0)
#     group_by(line_order) %>%
#     mutate(date_num = changeFun(as.numeric(date)) ) %>%
#     ungroup() %>%
#     # for every day calculate the amount to perturb a given line
#     group_by(date) %>%
#     mutate(pert_val = amp*sin(freq*date_num + offset) ) %>%
#     ungroup() %>%
#     # add perturbation to each line-day
#     mutate(pert_rv = rv_template_0.5 + pert_val) %>%
#     # include mean and median pert rv for each day
#     group_by(date) %>%
#     mutate(mean_pert_rv_day = mean(pert_rv),
#            med_pert_rv_day = median(pert_rv)) %>%
#     ungroup() %>%
#     # include mean and median pert rv for each line
#     group_by(line_order) %>%
#     mutate(mean_pert_rv_line = mean(pert_rv),
#            med_pert_rv_line = median(pert_rv)) %>%
#     ungroup()
# }

# ###########################
# ## create sparse matrix for the average over lines via left multiplying ##
# ###########################
# 
# # num_days : int for the number of days in the dataset
# # num_lines : int for the number of unique lineIDs
# 
# # return avg_matrix which is a num_days x ( num_lines*num_days ) sparse matrix
# # when used to left multiply a matrix with num_lines * num_days it returns the average over lines
# # assuming that this is grouped by individual days 
# avgMat = function(num_days, num_lines) {
#   # initialize an empty matrix
#   avg_matrix = Matrix(0, nrow = num_days, ncol = num_days * num_lines, sparse = T)
#   
#   # generate the average operator matrix
#   for (i in 1:num_days) {
#     avg_matrix[i, ((i-1) * num_lines + 1):(i * num_lines)] = 1/num_lines
#   }
#   
#   return(avg_matrix)
# }
# 
# ###########################
# ## get the empirical marginal mean for each day associated with each bootstrap coefficient ##
# ###########################
# 
# # bootcoef_df : R x p dataframe from casebootCoef, each row is a boot replicate and each column is a predictor
# # X_star : ( num_lines*num_days ) x p matrix
# # this is the "reference grid" where we have every line_order and date combination, and replace the continuous variables with their mean
# # ensure that the rows are ordered by date i.e the first 1:N values are for day 1, the N+1:2N are for day 2, and so on..
# # this ought to be a sparse design matrix for efficient computation 
# # num_lines and num_days are ints for the number of unique lines and days that we have in our dataframe
# 
# # returns day_empmean the matrix of num_days x R for the empirical marginal mean for the days
# empMeanDay = function(bootcoef_df, X_star, num_lines, num_days) {
#   # number of replicates and predictors
#   R = dim(bootcoef_df)[1]
#   p = dim(bootcoef_df)[2]
#   
#   if (dim(X_star)[2] != p) {
#     print("number of dimesnions of X_star doesn't match bootcoef_df")
#     return(NA)
#   } else if (!all( colnames(X_star) == colnames(bootcoef_df) )) {
#     print("column names do not match")
#     return(NA)    
#   }
#   
#   ## CHECK TO MAKE SURE WE'RE AVERAGING PROPERLY ##
#   
#   # calculate ( num_lines*num_days ) x R matrix of predictions associated with every row of the refGrid
#   pred_mat = X_star %*% t(bootcoef_df)
#   # calculate the matrix of num_days x R for the empirical marginal mean for the days
#   day_empmean = avgMat(num_days,num_lines) %*% pred_mat
#   
#   return(day_empmean)
# }
# # standardize our shape dataset
# 
# train_test_standardize = function(current_df, covars, test_index = c(1), response = "rv_template_0.5",
#                                   RESPONSE_CHANGE = T, TRAIN_TEST_SPLIT = T) {
#   
#   # modify response if we want to measure a change in response
#   if (RESPONSE_CHANGE) {
#     current_df = current_df %>%
#       group_by(line_order) %>%
#       mutate_at(response,changeFun) %>%
#       ungroup() %>%
#       rename_at(response, ~ paste0(., '_change'))
#   }
#   
#   # if we want to do training and testing, set TRAIN_TEST_SPLIT to T and this makes a train and test set,
#   # otherwise it just makes a train set of the current df
#   if (TRAIN_TEST_SPLIT) {
#     
#     # create training data without test_index
#     train = current_df[-test_index,]
#     # create testing data with test_index
#     test = current_df[test_index,]
#   } else {
#     train = current_df
#   }
# 
#   # get means of our covariates per line (for trained days)
#   train_means = train %>%
#     group_by(line_order) %>%
#     summarize_at(covars, mean)%>%
#     rename_at(covars, ~ paste0(., '_mean'))
#   # get sds of our centered covariates (for trained days)
#   train_sds = train %>%
#     group_by(line_order) %>%
#     mutate_at(covars,centerFun) %>%
#     ungroup() %>%
#     summarize_at(covars,sd)
#   
#   # for train data, center covars by line, scale by the entire column
#   # train = train %>%
#   #   group_by(line_order) %>%
#   #   mutate_at(covars,centerFun) %>%
#   #   ungroup() %>%
#   #   mutate_at(covars, scaleFun) %>%
#   #   rename_at(covars, ~ paste0(., '_center')) %>%
#   #   mutate(date = factor(date))
#   
#   train = train %>%
#     group_by(line_order) %>%
#     mutate_at(covars,centerFun) %>%
#     ungroup() %>%
#     mutate_at(covars, scaleFun) %>%
#     rename_at(covars, ~ paste0(., '_center'))
# 
#   # if we dont want to do training and testing, set TRAIN_TEST_SPLIT to F and this just returns the train set
#   if (!TRAIN_TEST_SPLIT) {
#     return(train)
#   }
#   
#   ## CHECK SEVERAL TIMES TO MAKE SURE THIS WORKS ##
#   # this section normalizes test data by means and sd of train data #
#   
#   # merge the training shape measurement means by line_order
#   test = merge(test, train_means, by = "line_order")
#   
#   # loop over covars, centering by the mean
#   for (covar in covars) {
#     test[[paste0(covar, "_center")]] = test[[covar]] - test[[paste0(covar, "_mean")]]
#   }
#   
#   # divide by standard deviation of train data
#   test[,str_c(covars,"_center")] = test[,str_c(covars,"_center")]/train_sds[rep(1,dim(test)[1]),]
#   
#   # drop covars and mean covars from test data
#   test = test[, !colnames(test) %in% c(covars, str_c(covars,"_mean") )]
#   
#   # set date as a factor
#   # test = test %>%
#   #   mutate(date = factor(date))
#   
#   return(list(train.df = train,
#               test.df = test))
# }
# ###########################
# ## case resampling for model coefficients see Algorithm 6.2 in Davison and Hinkley ##
# ###########################
# 
# # X (n x p) : model matrix
# # Y (n x 1) : responses
# # R : int for number of bootstrap samples to take
# 
# # returns an R x p matrix of the bootstrap coefficients
# casebootCoef2 = function(X, Y, R) {
#   START_TIME = Sys.time()
#   # initialize bootstrap coefficients as R x p matrix
#   coef_boot = matrix(nrow = R, ncol = dim(X)[2])
#   colnames(coef_boot) = colnames(X)
#   
#   # begin bootstrap
#   for (r in 1:R) {
#     # sample with replacement of indices 
#     idx_boot = sample( seq(1,nrow(X)) , replace = T)
#     # boot strap design matrix and response
#     X_boot = X[idx_boot,]
#     Y_boot = Y[idx_boot]
#     
#     XtX = t(X_boot) %*% X_boot
#     # check if our bootstrap design matrix is rank deficient, skip if so
#     if ( rankMatrix(XtX , method = "qr")[1] < dim(X)[2] ) {
#       print("warning: rank deficient bootstrap design matrix, skipping this sample")
#       next
#     }
#     
#     # get inverse of (X_boot^T X_boot)
#     XtX_inv = solve( XtX )
#     # get X^T Y
#     XtY = t(X_boot) %*% Y_boot
#     # store coefficient values
#     coef_boot[r,] = (XtX_inv %*% XtY)[,1]
#     
#     # print every 100th r
#     if (r%%100 == 0)  {
#       print(r)
#       print( Sys.time() - START_TIME )
#       }
# 
#     
#   }
#   
#   print( Sys.time() - START_TIME )
#   return(coef_boot)
#   
# }

# X (n x p) : model matrix
# Y (n x 1) : responses
# R : int for number of bootstrap samples to take

# returns an R x p matrix of the bootstrap coefficients
# casebootCoef = function(X, Y, R) {
#   START_TIME = Sys.time()
#   # initialize the indices used to generate bootstrap coefficients as R x n matrix
#   indices_boot = matrix(nrow = R, ncol = dim(X)[1])
#   # initialize bootstrap coefficients as R x p matrix
#   coef_boot = matrix(nrow = R, ncol = dim(X)[2])
#   colnames(coef_boot) = colnames(X)
#   # initialize bootstrap residuals as R x n matrix
#   resid_boot = matrix(nrow = R, ncol = dim(X)[1])
#   
#   # begin bootstrap
#   for (r in 1:R) {
#     # sample with replacement of indices 
#     idx_boot = sample( seq(1,nrow(X)) , replace = T)
#     # store indices
#     indices_boot[r,] = idx_boot
#     
#     # boot strap the design matrix and response
#     X_boot = X[idx_boot,]
#     Y_boot = Y[idx_boot]
#     
#     # fit the bootstrap matrices
#     fit_boot = sparseLM(X_boot, Y_boot, PRINT_TIME = F)
#     # store bootstrap coefficients
#     coef_boot[r,] = fit_boot$beta_hat[,1]
#     # store bootstrap residuals
#     resid_boot[r,] = fit_boot$resid[,1]
#     
#     # print every 100th r
#     if (r%%100 == 0)  {
#       print(r)
#       print( Sys.time() - START_TIME )
#     }
#     
#   }
#   
#   print( Sys.time() - START_TIME )
#   return( list(coef_boot = coef_boot,
#                resid_boot = resid_boot,
#                indices_boot = indices_boot) )
#   
# }
# ###########################
# ## Standardize our data, train-test split, and pca regression as options ##
# ###########################
# 
# # requires:
# # current_df - dataframe with lineID, date, and various covariates and a response as columns
# # covars - vector of strings for covariates that we want to standardize
# # covars_pca - optional vector of strings for the covariates that we want to do PCA on
# # test_index - vector of logical for the indices that we would like to split from training data
# # response - string for the response column
# # lineID - string for the lineID
# 
# standardizeWithPC = function(current_df, covariates = NULL, covars_pca = NULL, test_index = NULL, response = "rv_template_0.5", lineID = "line_order") {
#   
#   # if we want to do training and testing, include indices in test_index and this makes a train and test set,
#   # otherwise it just makes a train set of the current df
#   if (!is.null( test_index )) {
#     # create training data without test_index
#     train = current_df[!test_index,]
#     # create testing data with test_index
#     test = current_df[test_index,]
#   } else {
#     train = current_df
#   }
#   
#   # if no covariates are included, return original df
#   if (is.null(covariates)) {
#     return( list(train.df = train) )
#   }
#   
#   ## this section normalizes train data by means and sd of train data ##
#   
#   # get means of our covariates per line (for training days)
#   train_means = train %>%
#     group_by(!!sym(lineIDname)) %>%
#     summarize_at(covariates, mean)%>%
#     rename_at(covariates, ~ paste0(., '_mean'))
#   # get sds of our centered covariates (for training days)
#   train_sds = train %>%
#     group_by(!!sym(lineIDname)) %>%
#     mutate_at(covariates,centerFun) %>%
#     ungroup() %>%
#     summarize_at(covariates,sd)
#   # merge the training shape measurement means by lineIDname
#   train = merge(train, train_means, by = lineIDname)
#   # loop over covariates, centering by the mean
#   for (covar in covariates) {
#     train[[paste0(covar, "_centered")]] = train[[covar]] - train[[paste0(covar, "_mean")]]
#   }
#   # divide by standard deviation of train data
#   train[,str_c(covariates,"_centered")] = train[,str_c(covariates,"_centered")]/train_sds[rep(1,dim(train)[1]),]
#   # drop mean covariates from train data
#   train = train[, !colnames(train) %in% str_c(covariates,"_mean") ]
#   
#   ## this section does PCA on the select PCA covars with the train data ##
#   
#   if ( !is.null( covars_pca ) ) {
#     # select all covars to have PCA done on for train data
#     pca_covars_train = train %>%
#       select( all_of(str_c(covars_pca,"_centered")) )
#     # do pca on the train data (already centered and scaled)
#     pca_train = prcomp(pca_covars_train, retx = T)
#     # add PCs to training data
#     train = cbind(train, pca_train$x)
#   }
#   
#   ## if we dont want to do training and testing, this just returns the train set ##
#   if (is.null( test_index )) {
#     return(list(train.df = train,
#                 train_means = train_means,
#                 train_sds = train_sds))
#   }
# 
#   ## this section normalizes test data by means and sd of train data ##
#   
#   # merge the training shape measurement means by lineID
#   test = merge(test, train_means, by = lineID)
#   # loop over covars, centering by the mean
#   for (covar in covars) {
#     test[[paste0(covar, "_centered")]] = test[[covar]] - test[[paste0(covar, "_mean")]]
#   }
#   # divide by standard deviation of train data
#   test[,str_c(covars,"_centered")] = test[,str_c(covars,"_centered")]/train_sds[rep(1,dim(test)[1]),]
#   # drop mean covars from test data
#   test = test[, !colnames(test) %in% str_c(covars,"_mean") ]
#   
#   ## this section applies to test data the pca rotation##
#   
#   if ( !is.null( covars_pca ) ) {
#     # select all covars to have PCA done on for test data
#     pca_covars_test = test %>%
#       select( all_of(str_c(covars_pca,"_centered")) )
#     
#     # add PCs to test data using the training rotation
#     test = cbind(test, as.matrix(pca_covars_test) %*% pca_train$rotation )    
#   }
#   
#   return(list(train.df = train,
#               test.df = test,
#               train_means = train_means,
#               train_sds = train_sds))
# }