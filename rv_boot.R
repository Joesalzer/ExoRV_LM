source("readSpectra.R")
library(Matrix)
library(parallel)
library(pbmcapply)

# capture command line arguments
args = commandArgs(trailingOnly = TRUE)

# arguments passed to command line
# model name
MODEL_NAME = args[1]
NUM_BOOTS = args[2]

## variables below to be set by the user ##
# working directory with the data
WD_DATA = "/Users/josephsalzer/research/exostat/"
# names of the timeID, lineID, and timeGroupID
TIME_ID_NAME = "date"
LINE_ID_NAME = "line_order"

# directory containing the bootstrap files
bootstrap_dir = str_c(WD_DATA, "models/", MODEL_NAME, "/wild_bootstraps" )
# create directory called "bootstrap" in model name
if ( !dir.exists(bootstrap_dir) ) {
  dir.create(bootstrap_dir)
}

# read in the model fit
model_fit = readRDS(str_c(WD_DATA, "models/", MODEL_NAME, "/model.rds"))
# lm fit, df, design matrix and leverages
fit_lm = model_fit$fit_lm
rv_df = model_fit$df
designMat = model_fit$designMat
leverages = readRDS(str_c(WD_DATA, "models/", MODEL_NAME, "/leverages.rds" ))

# vec of LineIDs  in completeLines_df
lineIDs = rv_df %>% group_by(!!sym(LINE_ID_NAME)) %>% summarize(n. = n()) %>% pull(!!sym(LINE_ID_NAME))
# vec of timeIDs in completeLines_df
timeIDs = as.Date(rv_df %>% group_by(!!sym(TIME_ID_NAME)) %>% summarize(n. = n()) %>% pull(!!sym(TIME_ID_NAME)))

T_ = length(timeIDs)
L_ = length(lineIDs)

# get modified residuals
u_hat = fit_lm$resid/(1-leverages)

# assign blocks for the bootstrap (by line)
block_size = length(timeIDs)

rv_df = rv_df %>%
  arrange(!!sym(LINE_ID_NAME),!!sym(TIME_ID_NAME)) %>%
  mutate(boot_block = factor( ceiling(row_number() / block_size) ) )

wildBoot_parallel = function(seedID) {
  
  # set a given seed for replication
  set.seed(seedID)
  
  # ensure that the rv_df is arranged first by lineID and then timeID, and that these are factors
  rv_df = rv_df %>%
    arrange(!!sym(LINE_ID_NAME), !!sym(TIME_ID_NAME)) %>%
    mutate( lineID = factor(!!sym(LINE_ID_NAME)),
            timeID = factor(!!sym(TIME_ID_NAME)))
  
  # rademacher random variable for each block in the dataset
  v = sample( c(-1, 1), size = length(lineIDs), replace = TRUE)
  # generate wild bootstrap residuals by block
  u_star = u_hat*v[rv_df$boot_block]
  
  # construct bootstrap sample
  y_star = fit_lm$y_hat + u_star
  # re-estimate the model on bootstrap sample
  fitstar_lm = sparseLM(designMat, y_star, PRINT_TIME = F)
  
  # store data
  saveRDS(fitstar_lm$beta_hat[,1],
          file = str_c(bootstrap_dir, "/bootstraps_list_", seedID, ".rds" ) )
  
}

boot_straps <- pbmclapply(100*c(1:NUM_BOOTS), FUN = wildBoot_parallel, mc.cores = 7)