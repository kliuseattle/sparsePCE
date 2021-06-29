library(iterpc)
library(multipol)
library(groc)
library(glmnet)
library(MASS)
library(ggplot2)
library(dplyr)

poly_term <- function (p, d){
  I <- iterpc((p+1),d, ordered = T, replace = T) 
  temp <- getall(I)
  row <- rep(F, nrow(temp))
  for(i in 1:nrow(temp)){
    if(sum(temp[i,]-1) < (p+1) ){
      row[i] <- T
    }
  }
  res <- temp[row,]
  res <- res[,ncol(res):1]
  res <- res -1 
  #res <- apply(res, 1, rev)
  return(res)
}

create_X_feature <- function(X, p){
  d <- ncol(X)
  N <- nrow(X)
  polyorder_collect <- poly_term(p, d)[-1,]
  X_feature <- NA
  for(i in 1:nrow(polyorder_collect)){
    X_temp <- rep(1,N)
    for(j in 1:ncol(polyorder_collect)){
      X_temp <- X_temp * X[,j]^polyorder_collect[i,j]
    }
    X_feature <- cbind(X_feature, X_temp)
  }
  return(X_feature[,-1])
}

### comments: normalize the data first then select only feature column then transform.  
sparse_GS <- function(X, Y, thres_1){
  n <- nrow(X)
  m <- ncol(X)
  ortho_vec <- rep(1, n)/(n-1)
  temp_vec <- rep(1,n)/n
  flag <- TRUE
  select_index <- rep(0,m)
  pick_order <- 1 
  #X <- scale(X)
  while(flag == TRUE){
    #print(pick_order)
    X <- transform(X, temp_vec, select_index)
    
    
    select_res <- select(X, Y, thres_1, ortho_vec, select_index, pick_order)
    X <- select_res[[1]]
    temp_vec <- select_res[[2]]
    flag <- select_res[[3]]
    ortho_vec <- select_res[[4]]
    select_index <- select_res[[5]]
    pick_order <- pick_order + 1 
    

  }
  return(list('ortho_vec' = ortho_vec, 'select_idx' = select_index))
}

transform <- function(X, temp_vec, select_index){
  m <- ncol(X)
  n <- nrow(X)
  for(i in 1:m){
    if(select_index[i] ==0) {
      v <- X[,i]
      v <- v - c((X[,i] %*% temp_vec)/(temp_vec %*% temp_vec)) * temp_vec
      r_jj <- sqrt(sum(v^2))/sqrt((n-1))
      X[,i] <- v/r_jj
    }
  }
  return(X)
}



select <- function(X, y, thres_1, ortho_vec, select_index, pick_order){
  unselect <- which(select_index ==0 )
  corr_vec <- rep(0, ncol(X))
  if(length(unselect) ==0){
    flag = FALSE
    temp_vec = NA
    return(list(X, temp_vec, flag, ortho_vec, select_index))
  }
  else{
    for(i in 1:length(unselect)){
      corr_vec[unselect[i]] <- cor(X[,unselect[i]], y)
      if(abs(corr_vec[unselect[i]]) < thres_1){
        select_index[unselect[i]] <- -1 
      }
    }
    pick_one <- which(abs(corr_vec) == max(abs(corr_vec)))[1]
    #print(corr_vec)
    #print(pick_one)
    if(abs(corr_vec[pick_one]) >= thres_1){
      flag = TRUE
      temp_vec = X[,pick_one] 
      select_index[pick_one] <- pick_order
      ortho_vec <- cbind(ortho_vec, temp_vec)
      return(list(X, temp_vec, flag, ortho_vec, select_index))
    }  
    if(abs(corr_vec[pick_one]) < thres_1){
      flag = FALSE
      temp_vec = NA
      select_index[pick_one] <- -1
      return(list(X, temp_vec, flag, ortho_vec, select_index))
    }    
  }
}




