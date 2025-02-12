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

  ## if we dont want to do training and testing split, this just returns the train set ##
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
  
  # sum of squared residuals
  SSR = drop(crossprod(resid))
  
  # sigma2 hat and variance of beta_hat
  sigma2_hat = SSR / n 
  var_beta_hat = sigma2_hat * XtX_inv

  # total sum of squares
  SST = sum( ( Y - mean(Y) )^2 )
  # R-squared, AIC, BIC calculations
  R2 = 1-SSR/SST
  adj_R2 = 1 - (1- R2)*(n - 1)/(n - p - 1)
  AIC = n * log(SSR/n) + 2 * p
  BIC = n * log(SSR/n) + p * log(n)
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