require(tidyverse)
require(stringr)
require(dtw)

###########################
## function to save a distance between two rows of a matrix ##
###########################
# rows : a size 2 vector of which rows of the mat to take distance
# mat : a n x p matrix of n objects of size p that we want to take the distance of
# measure : string indicating what distance we want to calculate between rows, currently supporting euc, dtw, and ccf
dist_parallel = function(rows,mat,measure,...) {
  y1 = mat[rows[1],]
  y2 = mat[rows[2],]
  name1 = rownames(mat)[rows[1]]
  name2 = rownames(mat)[rows[2]]
  # create directory called "distances" in working directory
  if (!dir.exists( str_c(measure,"_distances") )) {dir.create(str_c(measure,"_distances"))}
  if (measure == "euc") {
    saveRDS(
      sqrt( sum( (y1-y2)^2 ) ),
      str_c(measure,"_distances/",name1,"#",name2,".rds")
    )
  } else if (measure == "ccf") {
    saveRDS(
      efficient_cross_correlation(y1,y2,...)$distance,
      str_c(measure,"_distances/",name1,"#",name2,".rds")
    )
  } else if (measure == "dtw") {
    saveRDS(
      dtw(y1,y2,...)$normalizedDistance,
      str_c(measure,"_distances/",name1,"#",name2,".rds")
    )
  }
}

###########################
## function that creates a list of all pairwise integers from 1 to n ##
###########################
# n - integer to create pairs from
create_pairs = function(n) {
  pairs_list = list()
  index = 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      pairs_list[[index]] = c(i, j)
      index = index + 1
    }
  }
  return(pairs_list)
}

###########################
## function to save a distance matrix from all pairwise distances ##
###########################
# mat : a n x p matrix of n objects of size p that we want to take the distance of
# measure : string indicating what distance we want to calculate between rows, currently supporting euc, dtw, and ccf
create_dist_mat = function(measure, mat) {
  # list of all pairwise distance files
  fileList = list.files(str_c(measure,"_distances"))
  # initialize matrices with 0's
  distance_matrix = matrix(0, nrow = nrow(mat), ncol = nrow(mat))
  rownames(distance_matrix) = rownames(mat)
  colnames(distance_matrix) = rownames(mat)
  
  # populate the matrix
  for (file in fileList) {
    # load the RDS file
    distance = readRDS( str_c(measure,"_distances/",file) )
    
    # extract days from the file name
    pair_vec = str_split_1(str_sub(file,end = -5),"#")
    
    # fill in the distance matrix (symmetrically)
    distance_matrix[pair_vec[1], pair_vec[2]] = distance
    distance_matrix[pair_vec[2], pair_vec[1]] = distance
    
    rm(distance,pair_vec,file)
  }
  
  saveRDS(distance_matrix, str_c(measure,"_distMat.rds"))
}

###########################
## euclidean distance between two vectors ##
###########################
euc_dist = function(y1, y2) sqrt(sum((y1 - y2) ^ 2))

###########################
## shift a logical vector by a specified lag ##
###########################
# vec : (logical) vector with a number of Falses followed by Trues followed by Falses
# lag : integer to shift the Trues by
shift_logical = function(vec, lag) {
  n = length(vec)
  if (lag > 0) {
    return( c( rep(F,lag), vec[1:(n-lag)] ) )
  } else if (lag < 0) {
    return( c( vec[(-lag + 1):n], rep(F,-lag) ) )
  } else {
    return(vec)
  }
}

###########################
## get cross correlation between two curves ##
###########################
# y1, y2 : vectors to take ccf of
# x_grid : a logical vector of which indices of y1, y2 to consider, should be a number of Falses followed by Trues followed by Falses
# max_lag : positive integer to get the dist between -max_lag, max_lag
cross_correlation = function(y1, y2, x_grid, max_lag) {
  lags = seq(-max_lag, max_lag)
  ccf_value = numeric(length(lags))
  for (i in seq_along(lags)) {
    lag = lags[i]
    ccf_value[i] = euc_dist(y1[x_grid], y2[shift_logical(x_grid,lag)])
  }
  return(data.frame(lag = lags, ccf_val = ccf_value))
}

###########################
## efficient ccf assuming the distance is convex, Eucledian dist ##
###########################
# y1, y2 : vectors to take ccf of
# x_grid : a logical vector of which indices of y1, y2 to consider, should be a number of Falses followed by Trues followed by Falses
# note that this function assumes that the ccf function is concave!!
efficient_cross_correlation = function(y1, y2, x_grid) {
  lag = 1
  left_ccf = euc_dist(y1[x_grid], y2[shift_logical(x_grid,-lag)])
  right_ccf = euc_dist(y1[x_grid], y2[shift_logical(x_grid,lag)])
  if (left_ccf < right_ccf) {
    while(left_ccf < right_ccf) {
      lag = lag + 1
      right_ccf = left_ccf
      left_ccf = euc_dist(y1[x_grid], y2[shift_logical(x_grid,-lag)])
    }
    return( list(shift = -lag+1,
                 distance = right_ccf) )
  } else {
    while(left_ccf > right_ccf) {
      lag = lag + 1
      left_ccf = right_ccf
      right_ccf = euc_dist(y1[x_grid], y2[shift_logical(x_grid,lag)])
    }
    return( list(shift = lag-1,
                 distance = left_ccf) )
  }
}

###########################
## Create omni matrix from n-by-n-by-m similarity array  ##
###########################
# M = omni(A)
# A : n-by-n-by-m array, A[,,i] is the adjacency matrix of i-th network.
# Returns
# M : nm-by-nm omnibus matrix
omni = function( A ) {
  # get size of array
  size_array = dim(A)
  # number of vertices
  nvx = size_array[2]
  # number of subjects(?)
  nsubj = size_array[3]
  
  # initialize omnibus matrix
  M = matrix(data=0,nrow=nsubj*nvx,ncol=nsubj*nvx)
  
  for (s in 1:nsubj) {
    for (t in s:nsubj) {
      i = nvx*(s-1) + 1
      j = nvx*(t-1) + 1
      B = (A[,,s] + A[,,t])/2
      M[ i:(i+nvx-1), j:(j+nvx-1) ] = B
      M[ j:(j+nvx-1), i:(i+nvx-1) ] = B
    }
  }
  
  # return error if M is not symmetric
  if (!isSymmetric(M)) {
    stop('Why did M not turn out symmetric?')
  }
  
  return(M)
  
}

###########################
## Create omnibus embedding from n-by-n-by-m array representing m n-by-n similarity matrices ##
###########################
# [X,evals] = OMNI_EMBED(A,d)
# Construct the d-dimensional omnibus embedding from array A.
# A : n-by-n-by-m array representing m n-by-n matrices.
# d : positive integer representing embedding dimension.
# Returns
# X : n-by-d matrix, i-th row is embedding of i-th vertex.
# evals : d-vector of eigenvalues of A, sorted descending.
omni_embed = function( A, d ) {
  # get omnibus matrix from A
  M = omni(A)
  # now do eigendecomp of M
  decompM = eigen(M)
  # diagonal matrix containing the d largest eigenvalues on the main diagonal
  S = diag(decompM$values[1:d])
  # matrix V whose columns are the corresponding eigenvectors
  V = decompM$vectors[,1:d]
  # embedding
  X = V %*% sqrt(S)
  
  return( list("embedding" = X, "evals" = decompM$values[1:d]) )
  
}

###########################
## Create adjecency embedding from n-by-n-by-m array representing an n-by-n similarity matrix ##
###########################
# Construct the d-dimensional adj embedding from matrix M.
# M : n-by-n array
# d : positive integer representing embedding dimension.
# Returns
# X : n-by-d matrix, i-th row is embedding of i-th vertex.
# evals : d-vector of eigenvalues of A, sorted descending.
adj_embed = function( M, d ) {
  # do eigendecomp of M
  decompM = eigen(M)
  # diagonal matrix containing the d largest eigenvalues on the main diagonal
  S = diag(decompM$values[1:d])
  # matrix V whose columns are the corresponding eigenvectors
  V = decompM$vectors[,1:d]
  # embedding
  X = V %*% sqrt(S)
  
  return( list("embedding" = X, "evals" = decompM$values[1:d]) )
  
}

###########################
## Two sample T2 test ##
###########################
# embedding : a dataframe with a network_id column indicating the two networks, and columns X1, X2, ..., X(ndim)
TwoSampleT2Test = function(embedding, ndim = 3) {

  # get each network label
  network_ids = unique( embedding$network_id )
  # get each n (number of lambda-orders) x d (dimension embedding) dataframes
  X = embedding %>% filter(network_id == network_ids[1]) %>% select( str_c("X",c(1:ndim)) )
  Y = embedding %>% filter(network_id == network_ids[2]) %>% select( str_c("X",c(1:ndim)) )

  nx = nrow(X)
  ny = nrow(Y)
  delta = colMeans(X) - colMeans(Y)
  p = ncol(X)
  Sx = cov(X)
  Sy = cov(Y)
  S_pooled = ((nx-1)*Sx + (ny-1)*Sy)/(nx+ny-2)
  t_squared = (nx*ny)/(nx+ny) * t(delta) %*% solve(S_pooled) %*% (delta)
  statistic = t_squared * (nx+ny-p-1)/(p*(nx+ny-2))
  p_value = pf(statistic, p, nx+ny-p-1, lower.tail=FALSE)
  #print(glue("Test statistic: {statistic}
  # Degrees of freedom: {p} and {nx+ny-p-1}
  # p-value: {p_value}"))
  return(list(TestStatistic=statistic, p_value=p_value))
}
# ## build a distance matrix, each row/col will be a pair of line_order-variable##
# 
# # vars : d-vector of strings, vars[i] is the ith var
# # select_df : dataframe of observations to build the corr array, n line_orders m dates and contains variables vars
# # measure: string of what dist measurement to use between time series, default is pearson correlation
# 
# # Returns
# # dist_matrix : (n x num_vars) x (n x num_vars) matrix where each row/col is a line_order-variable pair
# 
# build_dist_matrix = function( select_df, vars, lineID = "line_order", timeID = "date", measure = "euclidean" ) {
#   
#   # start the clock
#   ptm = proc.time()
# 
#   # initialize list of dataframes 
#   vars_list = list()
#   
#   # pivot the dataframe for each variable, loop over all variables
#   for (i in 1:length(vars)) {
#     cat("variable number: ", i, "\n")
#     var = vars[i]
#     
#     vars_list[[i]] = select_df %>%
#       dplyr::select(line_order, date, sym(var)) %>%
#       mutate( line_order = str_c(var,":",line_order) ) %>%
#       pivot_wider(names_from = line_order, values_from = sym(var)) %>%
#       select(!date)
#     cat("\t time after pivot wider: ", (proc.time() - ptm)[1], "\n")
#   }
#   
#   # bind togehter vars list by column
#   vars_df = do.call(cbind, vars_list)
#   rm(vars_list, var, i)
#   cat("bind together list of var dfs: ", (proc.time() - ptm)[1], "\n")
#   
#   if (measure == "dtw") {
#     dist_matrix = t( vars_df %>% scale() ) %>% dist(method = "dtw") %>% as.matrix()
#   } else if (measure == "euclidean") {
#     dist_matrix = t( vars_df %>% scale() ) %>% dist() %>% as.matrix()
#   }
#   
#   cat("time building distance matrix: ", (proc.time() - ptm)[1], "\n")
#   
#   return(dist_matrix)
#   
# }
# ## convert our large similarity matrix into the arrays of the individual correlation matrices for single variables##
# 
# # num_vars : integer of number of variables under consideration
# # mat : (num_lines x num_vars) x (num_lines x num_vars) matrix where each row/col is a line_order-variable pair
# 
# # Returns
# # multi_array : num_lines x num_lines x num_vars array where each obs is a correlation matrix
# 
# mat_to_array = function(mat, num_vars) {
#   # number of lines in the matrix
#   num_lines = dim(mat)[1]/num_vars
#   # initialize mutli-dimensional array of correlations
#   multi_array = array(NA,c(num_lines,num_lines,num_vars))
#   # vector of line names
#   line.names = unique( str_split_fixed(rownames(mat), ":", n = 2)[,2] )
#   
#   for (i in 1:num_vars) {
#     multi_array[,,i] = mat[ (num_lines*(i-1)+1):(i*num_lines), (num_lines*(i-1)+1):(i*num_lines) ]
#     
#     # name rows and cols
#     if (i == num_vars) {
#       colnames( multi_array ) = line.names
#       rownames( multi_array ) = line.names
#     }
#     
#   }
#   
#   return(multi_array)
# }
# ## function to build similarity array ##
# # vars : d-vector of strings, vars[i] is the ith var
# # select.df : dataframe of observations to build the corr array, n line_orders m dates 
# # measure: string of what similarity measurement to use between time series, default is pearson correlation
# # Returns
# # multi_network : n x n x d array where each obs is a correlation matrix
# build_sim_array = function( select.df, vars, measure = "pearson" ) {
# 
#   # number of variables
#   n_vars = length(vars)
#   # number of line_order pairs
#   n_lams = dim(select.df %>% count(line_order))[1]
#   # initialize mutli-dimensional array of correlations
#   multi_network = array(NA,c(n_lams,n_lams,n_vars))
#   
#   # start the clock
#   ptm = proc.time()
#   
#   for (i in 1:n_vars) {
#     cat("variable number: ", i, "\n")
#     
#     # select variable
#     var = vars[i]
#     
#     # pivot original dataframe to wide to get variable measurment over time
#     var.df = select.df %>%
#       dplyr::select(line_order, date, sym(var)) %>%
#       pivot_wider(names_from = line_order, values_from = sym(var)) %>%
#       select(!date)
#     
#     cat("\t time after pivot wider: ", (proc.time() - ptm)[1], "\n")
#     
#     # put correlation matrix into our multidim network
#     
#     if (measure == "pearson") {
#       multi_network[,,i] = abs( cor( var.df  ) )
#     } else if (measure == "spearman") {
#       multi_network[,,i] = abs( cor( var.df, method = "spearman"  ) )
#     } else if (measure == "dtw") {
#       multi_network[,,i] = t( var.df %>% scale() ) %>% dist(method = "dtw") %>% as.matrix() %>% pr_dist2simil() 
#     } else if (measure == "euclidean") {
#       multi_network[,,i] = t( var.df %>% scale() ) %>% dist() %>% as.matrix() %>% pr_dist2simil() 
#     } else if (measure == "dtw_noscale") {
#       multi_network[,,i] = t( var.df ) %>% dist(method = "dtw") %>% as.matrix() %>% pr_dist2simil() 
#     }
#     
#     cat("\t time building similarity matrix: ", (proc.time() - ptm)[1], "\n")
#     
#     # name rows and cols
#     if (i == n_vars) {
#       colnames( multi_network ) = colnames( var.df )
#       rownames( multi_network ) = colnames( var.df )
#     }
#     
#     rm(var.df)
#     
#   }
#   
#   return(multi_network)
#   
# }
# ## permute a similarity matrix, each row/col will be a pair of line_order-variable##
# 
# # sim.mat : (n x 2) x (n x 2) matrix where each row/col is a line_order-variable pair, only works for 2 variables
# # indices : vector of lines (from 1 to n) to permute
# 
# # Returns
# # permuted.mat : the permuted similarity matrix
# 
# permute_sim_matrix = function(sim.mat, indices) {
#   # number of lines
#   n = dim(sim.mat)[1]/2
#   # initialize permuted matrix
#   permuted.mat = sim.mat
#   for (idx in indices) {
#     # Swap the rows + columns corresponding to the same line
#     permuted.mat[, c(idx, idx + n)] = permuted.mat[, c( idx + n, idx )]
#     permuted.mat[c(idx, idx + n), ] = permuted.mat[c( idx + n, idx ), ]
#   }
#   
#   # return error if permuted.mat is not symmetric
#   if (!isSymmetric(permuted.mat)) {
#     stop('Why did the permuted matrix not turn out symmetric?')
#   } else if ( all(diag(permuted.mat) != 1) ) {
#     stop('Why does the permuted matrix not have ones on the diagonals?')
#   }
#   
#   return(permuted.mat)
#   
# }
# ## function to build correlation array ##
# 
# # vars : d-vector of strings, vars[i] is the ith var
# # select.df : dataframe of observations to build the corr array, n line_orders m dates 
# # Returns
# # multi_network : n x n x d array where each obs is a correlation matrix
# build_cor_array = function( vars, select.df ) {
#   # number of variables
#   n_vars = length(vars)
#   # number of line_order pairs
#   n_lams = dim(select.df %>% count(line_order))[1]
#   # initialize mutli-dimensional array of correlations
#   multi_network = array(NA,c(n_lams,n_lams,n_vars))
#   
#   for (i in 1:n_vars) {
#     # select variable
#     var = vars[i]
#     
#     # pivot original dataframe to wide to get variable measurment over time
#     var.df = select.df %>%
#       dplyr::select(line_order, date, sym(var)) %>%
#       pivot_wider(names_from = line_order, values_from = sym(var)) %>%
#       select(!date)
#     
#     # put correlation matrix into our multidim network
#     
#     multi_network[,,i] = abs( cor( var.df  ) )
#     #multi_network[,,i] = cor( var.df  )
#     #multi_network[,,i] = abs( cor( var.df, method = "spearman"  ) )
#     
#     # name rows and cols
#     if (i == n_vars) {
#       colnames( multi_network ) = colnames( var.df )
#       rownames( multi_network ) = colnames( var.df )
#     }
#     
#     rm(var.df)
#     
#   }
#   
#   return(multi_network)
#   
# }
# ## Permute dataframe ##
# 
# # vars : d-vector of strings, vars[i] is the ith var
# # select.df : dataframe of observations to build the corr array, n line_orders m dates 
# # Returns
# # new dataframe where permvar1 and permvar2 are permuted versions of the vars 
# permuteDF = function(vars, select.df) {
#   
#   # modify dataframe for whether we will permute our column variables for a vertex/line_order for ~every date~
#   perm.df = select.df %>%
#     dplyr::select(line_order, date, all_of(vars)) %>%
#     group_by(line_order) %>%
#     mutate(toPerm = rbinom(n = 1, size = 1, prob = .5) ) %>%
#     ungroup()
#   
#   # permute columns if toPerm==1, otherwise keep columns the same
#   return( perm.df %>%
#             mutate( permvar1 = ifelse(toPerm==1, !!sym(vars[2]), !!sym(vars[1])  ),
#                     permvar2 = ifelse(toPerm==1, !!sym(vars[1]), !!sym(vars[2])  ) )
#   )
# }
# 
# ## Create distance matrix between embeddings from omni or adjacency embedding ##
# 
# # embedding: dataframe of each line_order-networkID pair, where there are columns for line_order, network_id, and each of the embedding dimensions labeled X1 through X(n_dim)
# # ndim: integer for the size of the dimension of embedding
# # measure: string for which distance measure to use, defaults to euclidean
# 
# # returns the distance matrix with rows and columns labeled by line_order-networkID
# 
# dist_matrix = function(embedding, n_dim, measure = "euclidean") {
#   # create distance matrix
#   dist.matrix = embedding %>%
#     select(str_c("X",seq(1,n_dim))) %>%
#     stats::dist(method = measure) %>%
#     as.matrix( )
#   rownames(dist.matrix) = embedding %>% mutate(id = str_c(line_order,"-",network_id)) %>% pull(id)
#   colnames(dist.matrix) = rownames(dist.matrix)
#   
#   return(dist.matrix)
# }