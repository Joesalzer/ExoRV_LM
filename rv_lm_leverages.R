source("readSpectra.R")
library(Matrix)
library(parallel)
library(pbmcapply)

# capture command line arguments
args = commandArgs(trailingOnly = TRUE)

# model name
MODEL_NAME = args[1]
BATCH_SIZES = as.numeric(args[2])

## USER DEFINED VARS ##
# working directory with the data
WD_DATA = "/Users/josephsalzer/research/exostat/"
# read in the model fit
model_fit = readRDS(str_c(WD_DATA, "models/", MODEL_NAME, "/model.rds"))
## USER DEFINED VARS ##

# fit lm and design matrices
fit_lm = model_fit$fit_lm
X = model_fit$designMat
XtX_inv = (fit_lm$var_beta_hat)/fit_lm$sigma2_hat

# directory containing the leverage files
leverage_dir = str_c(WD_DATA, "models/", MODEL_NAME, "/leverages")
# create directory called "leverages" in model name
if ( !dir.exists(leverage_dir) ) {
  dir.create(leverage_dir)
}

# optimized leverage computation
compute_leverage_batches = function(indices) {
  start_row = indices[1]
  end_row = indices[2]
  
  # extract rows in a batch
  rows = X[start_row:end_row, , drop = FALSE]
  # compute leverage for all rows in the batch
  leverage = rowSums(rows * (rows %*% XtX_inv))
  
  names(leverage) = str_c("row",start_row:end_row)
  
  # save results in a single file for the batch
  batch_filename = file.path(leverage_dir, str_c("leverage_", start_row, "_", end_row, ".rds"))
  saveRDS(leverage, file = batch_filename)
}

# function to create batches
create_batches = function(n_rows, batch_size) {
  # start indices of batches
  starts = seq(1, n_rows, by = batch_size)
  # end indices of batches
  ends = pmin(starts + batch_size - 1, n_rows)
  # combine into a list of start and end indices
  batches = Map(c, starts, ends)
  return(batches)
}

# create batch indices
batches = create_batches(nrow(X), BATCH_SIZES)

# parallel computation of leverage
result <- pbmclapply(batches, FUN = compute_leverage_batches, mc.cores = 7)

## code below adds leverages to the model.rds, but produced some errors ##

# # get the list of all .rds files in the leverage directory
# leverage_files = list.files(path = leverage_dir, pattern = "\\.rds$", full.names = TRUE)
# # initialize an empty vector to store the combined leverage values
# leverages = vector()
# for (file in leverage_files) {
#   leverages = c(leverages, readRDS(file))
# }
# rm(leverage_files)
# 
# # add leverages to model_fit
# model_fit$leverages = leverages
# saveRDS(model_fit, file = str_c(WD_DATA, "models/", MODEL_NAME, "/model.rds" ) )