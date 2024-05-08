
###########################
# helper functions, to read in NEID datafiles and fit/evaluate models
###########################

require(tidyverse)
require(stringr)
require(rhdf5)

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
#

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
## Standardize our data, train-test split, and pca as options ##
###########################

# requires:
# current_df - dataframe of line_order, date, and various covariates and responses
# covars - vector of strings for covariates that we want to standardize (pass more than 1)
# covars_pca - optional vector of strings for the covariates that we want to do PCA on
# test_index - vector of ints for the indices that we would like to split from training data
# response - string for the response column
# RESPONSE_CHANGE - logical value for whether we want to mutate the response to change from the first day
# TRAIN_TEST_SPLIT - logical value for whether we want to do train test splitting or not

standardize = function(current_df, covars, covars_pca = NULL, test_index = NULL, response = "rv_template_0.5",
                                  RESPONSE_CHANGE = F, TRAIN_TEST_SPLIT = F) {
  # modify response if we want to measure a change in response
  if (RESPONSE_CHANGE) {
    current_df = current_df %>%
      group_by(line_order) %>%
      mutate_at(response,changeFun) %>%
      ungroup() %>%
      rename_at(response, ~ paste0(., '_change'))
  }
  
  # if we want to do training and testing, set TRAIN_TEST_SPLIT to T and this makes a train and test set,
  # otherwise it just makes a train set of the current df
  if (TRAIN_TEST_SPLIT) {
    # create training data without test_index
    train = current_df[-test_index,]
    # create testing data with test_index
    test = current_df[test_index,]
  } else {
    train = current_df
  }
  
  # get means of our covariates per line (for trained days)
  train_means = train %>%
    group_by(line_order) %>%
    summarize_at(covars, mean)%>%
    rename_at(covars, ~ paste0(., '_mean'))

  # get sds of our centered covariates (for trained days)
  train_sds = train %>%
    group_by(line_order) %>%
    mutate_at(covars,centerFun) %>%
    ungroup() %>%
    summarize_at(covars,sd)

  # modify training data to center/scale and to include original data
  train = train %>%
    group_by(line_order) %>%
    mutate_at(vars( all_of(covars) ), list(centered = ~ centerFun(.)) ) %>%
    ungroup() %>%
    mutate_at(vars(all_of(str_c(covars,"_centered"))), scaleFun)

  ## this section does PCA on the select PCA covars with the train data ##
  if ( !is.null( covars_pca ) ) {
    # select all covars to have PCA done on for train data
    
    pca_covars_train = train %>%
      select( all_of(str_c(covars_pca,"_centered")) )
    
    # do pca on the train data (already centered and scaled)
    pca_train = prcomp(pca_covars_train, retx = T)
    
    # add PCs to training data
    train = cbind(train, pca_train$x)
  }
  
  ## if we dont want to do training and testing, set TRAIN_TEST_SPLIT to F and this just returns the train set ##
  if (!TRAIN_TEST_SPLIT) {
    return(train)
  }
  
  ### CHECK SEVERAL TIMES TO MAKE SURE THIS WORKS ###
  
  ## this section normalizes test data by means and sd of train data ##
  
  # merge the training shape measurement means by line_order
  test = merge(test, train_means, by = "line_order")
  
  # loop over covars, centering by the mean
  for (covar in covars) {
    test[[paste0(covar, "_centered")]] = test[[covar]] - test[[paste0(covar, "_mean")]]
  }
  
  # divide by standard deviation of train data
  test[,str_c(covars,"_centered")] = test[,str_c(covars,"_centered")]/train_sds[rep(1,dim(test)[1]),]
  
  # drop mean covars from test data
  test = test[, !colnames(test) %in% str_c(covars,"_mean") ]
  
  ## this section applies to test data the pca rotation##
  
  if ( !is.null( covars_pca ) ) {
    # select all covars to have PCA done on for test data
    pca_covars_test = test %>%
      select( all_of(str_c(covars_pca,"_centered")) )
    
    # add PCs to test data using the training rotation
    test = cbind(test, as.matrix(pca_covars_test) %*% pca_train$rotation )    
  }
  
  return(list(train.df = train,
              test.df = test))
}

###########################
## fit a LM using model matrix and responses##
###########################

# X (n x p) : model matrix
# Y (n x 1) : responses

# returns a list of model fits

sparseLM = function(X, Y) {
  START_TIME = Sys.time()
  
  # get (X^T X)^-1 matrix
  XtX_inv = solve( t(X) %*% X )
  
  # get X^T Y
  XtY = t(X) %*% Y

  # beta hat, linear model coefficients
  beta_hat = XtX_inv %*% XtY

  # fitted value
  y_hat = X%*%beta_hat

  # residuals
  resid = Y - y_hat

  # number of dimensions
  p = dim(X)[2]

  # sigma2 hat, 
  sigma2_hat = t(resid) %*% resid/(length(Y) - p)

  # covariance matrix of coefficients
  var_beta_hat = drop(sigma2_hat)*XtX_inv
  
  print( Sys.time() - START_TIME )
  
  return( list( beta_hat = beta_hat,
                y_hat = y_hat,
                resid = resid,
                sigma2_hat = sigma2_hat,
                var_beta_hat = var_beta_hat))
}

###########################
## case resampling for model coefficients see Algorithm 6.2 in Davison and Hinkley ##
###########################

# X (n x p) : model matrix
# Y (n x 1) : responses
# R : int for number of bootstrap samples to take

# returns an R x p matrix of the bootstrap coefficients
casebootCoef = function(X, Y, R) {
  START_TIME = Sys.time()
  # initialize bootstrap coefficients as R x p matrix
  coef_boot = matrix(nrow = R, ncol = dim(X)[2])
  colnames(coef_boot) = colnames(X)
  
  # begin bootstrap
  for (r in 1:R) {
    # sample with replacement of indices 
    idx_boot = sample( seq(1,nrow(X)) , replace = T)
    # boot strap design matrix and response
    X_boot = X[idx_boot,]
    Y_boot = Y[idx_boot]
    
    XtX = t(X_boot) %*% X_boot
    # check if our bootstrap design matrix is rank deficient, skip if so
    if ( rankMatrix(XtX , method = "qr")[1] < dim(X)[2] ) {
      print("warning: rank deficient bootstrap design matrix, skipping this sample")
      next
    }
    
    # get inverse of (X_boot^T X_boot)
    XtX_inv = solve( XtX )
    # get X^T Y
    XtY = t(X_boot) %*% Y_boot
    # store coefficient values
    coef_boot[r,] = (XtX_inv %*% XtY)[,1]
    
    # print every 100th r
    if (r%%100 == 0)  {
      print(r)
      print( Sys.time() - START_TIME )
      }

    
  }
  
  print( Sys.time() - START_TIME )
  return(coef_boot)
  
}

###########################
## create sparse matrix for the average over lines via left multiplying ##
###########################

# num_days : int for the number of days in the dataset
# num_lines : int for the number of unique lineIDs

# return avg_matrix which is a num_days x ( num_lines*num_days ) sparse matrix
# when used to left multiply a matrix with num_lines*num_days it returns the average over lines
# assuming that this is grouped by individual days 
avgMat = function(num_days, num_lines) {
  # initialize an empty matrix
  avg_matrix = Matrix(0, nrow = num_days, ncol = num_days * num_lines, sparse = T)
  
  # generate the average operator matrix
  for (i in 1:num_days) {
    avg_matrix[i, ((i-1) * num_lines + 1):(i * num_lines)] = 1/num_lines
  }
  
  return(avg_matrix)
}

###########################
## get the empirical marginal mean for each day associated with each bootstrap coefficient ##
###########################

# bootcoef_df : R x p dataframe from casebootCoef, each row is a boot replicate and each column is a predictor
# X_star : ( num_lines*num_days ) x p matrix
# this is the "reference grid" where we have every line_order and date combination, and replace the continuous variables with their mean
# ensure that the rows are ordered by date i.e the first 1:N values are for day 1, the N+1:2N are for day 2, and so on..
# this ought to be a sparse design matrix for efficient computation 
# num_lines and num_days are ints for the number of unique lines and days that we have in our dataframe

# returns day_empmean the matrix of num_days x R for the empirical marginal mean for the days
empMeanDay = function(bootcoef_df, X_star, num_lines, num_days) {
  # number of replicates and predictors
  R = dim(bootcoef_df)[1]
  p = dim(bootcoef_df)[2]
  
  if (dim(X_star)[2] != p) {
    print("number of dimesnions of X_star doesn't match bootcoef_df")
    return(NA)
  } else if (!all( colnames(X_star) == colnames(bootcoef_df) )) {
    print("column names do not match")
    return(NA)    
  }
  
  ## CHECK TO MAKE SURE WE'RE AVERAGING PROPERLY ##
  
  # calculate ( num_lines*num_days ) x R matrix of predictions associated with every row of the refGrid
  pred_mat = X_star %*% t(bootcoef_df)
  # calculate the matrix of num_days x R for the empirical marginal mean for the days
  day_empmean = avgMat(num_days,num_lines) %*% pred_mat
  
  return(day_empmean)
}

###########################
## inject a planet of amplitude and frequency into our dataframe ##
###########################

# df : datframe that includes line_order, date, and rv_template_0.5
# amp : amplitude of the sin wave (defaults to 1)
# freq : frequency of rotation (defaults to earth-like rotation of 366 days)
# offset : horizontal shift of the planet signal (defaults to none)

injectPlanet = function(df, amp = 1, freq = 2*pi/366, offset = 0) {
  df %>%
    # for every line, get the date as a day value (starting at 0)
    group_by(line_order) %>%
    mutate(date_num = changeFun(as.numeric(date)) ) %>%
    ungroup() %>%
    # for every day calculate the amount to perturb a given line
    group_by(date) %>%
    mutate(pert_val = amp*sin(freq*date_num + offset) ) %>%
    ungroup() %>%
    # add perturbation to each line-day
    mutate(pert_rv = rv_template_0.5 + pert_val) %>%
    # include mean and median pert rv for each day
    group_by(date) %>%
    mutate(mean_pert_rv_day = mean(pert_rv),
           med_pert_rv_day = median(pert_rv)) %>%
    ungroup() %>%
    # include mean and median pert rv for each line
    group_by(line_order) %>%
    mutate(mean_pert_rv_line = mean(pert_rv),
           med_pert_rv_line = median(pert_rv)) %>%
    ungroup()
}



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
