library(GpGp)

permute_data = function(item,index){
  if (is.matrix(item)) {
    # If it's a matrix, subset the specified rows
    if(dim(item)[2]>1){
      item[index, , drop = FALSE] # Use drop=FALSE to prevent collapsing to vector if only one row/column remains  
    }
    else {
      # print('not matrix')
      item[index]
    }
    
  } else {
    # print('not matrix')
    # If it's NOT a matrix (e.g., a vector, data frame, etc.), keep it as is
    item[index]
  }
}


reorder_data = function(data, order = 'coordinate'){
  if(order == 'coordinate'){
    ## find the new order
    ord = order_coordinate(locs = as.matrix(data[,c('coord_x','coord_y')]))
    ## reorder all the elements in data
    ## lapply returns a list, we need to transform back to data.frame
    data = as.data.frame(lapply(data,permute_data,index = ord))
    return(data)
  } else if (order == 'maxmin'){
    ## find the new order
    ord = order_maxmin(locs = as.matrix(data[,c('coord_x','coord_y')]))
    ## reorder all the elements in data
    ## lapply returns a list, we need to transform back to data.frame
    data = as.data.frame(lapply(data,permute_data,index = ord))
    return(data)
  }
}

precision_column_calculation = function(col_index, Linv, NNarray){
  e = ifelse(1:dim(Linv)[1] == col_index, 1, 0)
  precision_column = Linv_t_mult(Linv = Linv,z = Linv_mult(Linv = Linv,z = e,NNarray = NNarray),NNarray = NNarray)
  return(precision_column)
}

# the data has to be ordered before 
# y_obs is the residuals of fitted outcome that is used to do kriging
leave_one_out_kriging = function(locs, y_obs, gp_model, gp_params, order = c("coordinate", "maxmin")){
  
  n = nrow(locs)
  locs = as.matrix(locs)
  
  ## re-order data
  order = match.arg(order)
  if(order == "coordinate"){
    ord = order_coordinate(locs = locs)  
  } else if(order == "maxmin"){
    ord = order_maxmin(locs = locs)
  }
  locs = permute_data(locs, ord)
  y_obs = permute_data(y_obs, ord)
  
  ## find the nearest neighbors
  NNarray = find_ordered_nn(locs = locs, m=30)
  ## calculate Linv
  Linv = vecchia_Linv(covparms = gp_params, covfun_name = gp_model, locs = locs, NNarray)
  
  # Then, calculate E[U_i|U_{-i}] for every i.
  y_pred = rep(NA,n)
  # should not include the position itself
  for(i in 1:n){
    if(i%%500==0)print(i)
    col_index = i
    precision_column = precision_column_calculation(col_index = col_index, Linv = Linv, NNarray = NNarray)
    y_pred[i] = -  sum(precision_column[ -col_index] * y_obs[- col_index ])/ precision_column[col_index]
  }
  
  # ## get full precision matrix
  # Precision = Linv %*% t(Linv)   # since Linv is triangular factor of precision
  # 
  # ## compute predictions in vectorized form
  # # diag_vec = diagonal elements P_ii
  # diag_vec = diag(Precision)
  # # numerator = (Precision %*% y) - P_ii * y_i
  # numer = Precision %*% y_obs - diag_vec * y_obs
  # y_pred = - numer / diag_vec
  
  
  # reverse to the original order
  inv_ord <- integer(length(ord))
  inv_ord[ord] <- 1:n
  y_pred <- permute_data(y_pred, inv_ord)
  
  return(y_pred)
}
