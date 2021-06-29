### Example 1 Ishigami function 
library(lhs)
library(lars)
data_create <- function(seed, size, LHS = FALSE){
  set.seed(seed)
  n <- size 
  if(LHS == FALSE){
    X1 <- runif(n, -pi, pi)
    X2 <- runif(n, -pi, pi)
    X3 <- runif(n, -pi, pi)
    Y <- sin(X1) + 7*sin(X2)^2 + 0.1*X3^4*sin(X1)
    return(list('X' = cbind(X1, X2, X3) , 'Y'=Y ))     
  }
 else{
   X_LHS <- randomLHS(n, 3)
   X1 <- qunif(X_LHS[,1], min= -pi, max= pi)
   X2 <- qunif(X_LHS[,2], min= -pi, max= pi)
   X3 <- qunif(X_LHS[,3], min= -pi, max= pi)
   Y <- sin(X1) + 7*sin(X2)^2 + 0.1*X3^4*sin(X1)
   return(list('X' = cbind(X1, X2, X3) , 'Y'=Y ))  
 }
}



mom_approxi_propose <- function(X, Y, p, thres){
  X_feature <- create_X_feature(X, p)
  res <- sparse_GS(X_feature, Y, thres)
  
  m1 <- lm(Y ~res$ortho_vec[,-1])
  coef <- m1$coefficients[!is.na(m1$coefficients)]
  mean_approx <- coef[1]
  sd_approx <- sqrt(sum(coef[-1]^2))
  return(list('mean_approx' = mean_approx, 'sd_approx'= sd_approx))
}



find_threshold <- function(X, Y, p, thres_list){
  n <- length(thres_list)
  valid_error <- rep(0, n)
  N <- nrow(X)
  d <- ncol(X)
  error_list <- rep(0, n)
  for(i in 1:n){
    #print(i)
    set.seed(123)
    thres <- thres_list[i]
    rep_order <- sample(N, N)
    temp_error <- 0
    for(j in 1:5){
      samp_j <- rep_order[((j-1)*N/5+1):(j*N/5)]
      X_train <- X[-samp_j, ]
      Y_train <- Y[-samp_j]
      X_test <- X[samp_j,]
      Y_test <- Y[samp_j]
      X_feature <- create_X_feature(X_train, p)
      res <- sparse_GS(X_feature, Y_train, thres)
      X_feature_full <- create_test_X_feature(X, res$select_idx, p)
      X_ortho_full <- gramschmidt(scale(X_feature_full))$Q
      XY_data <- as.data.frame(cbind(X_ortho_full, Y))
      X_ortho_train_data <- XY_data[-samp_j,]
      X_ortho_test_data <- XY_data[samp_j,]
      m1 <- lm(Y~., data = X_ortho_train_data)
      #print(summary(m1))
      y_hat <- predict(m1, X_ortho_test_data)
      temp_error <- temp_error + sum((Y_test - y_hat)^2)/length(Y_test)
    }
    error_list[i] <- temp_error/5
  }
  return(error_list)
}

create_test_X_feature <- function(X, select_idx, p){
  d <- ncol(X)
  poly_orders <- poly_term(p, d)[-1,]
  n <- max(select_idx)
  X_test <- NA
  N <- nrow(X)
  for (i in 1:n){
    #print(i)
    X_temp <- rep(1, N)
    temp_poly <- poly_orders[select_idx==i,]
    for(j in 1:length(temp_poly)){
      X_temp <- X_temp*X[,j]^temp_poly[j]
    }
    X_test <- cbind(X_test, X_temp)
  }
  X_test <- X_test[,-1]
  return(X_test)
}

### Find best threshold value from p = 1 to p = 9
data <- data_create(3, 200, LHS = TRUE)
X <- data$X
Y <- data$Y

var <- 49/8 + 0.1*pi^4/5 + 0.01*pi^8/18 + 0.5
sd_true <- sqrt(var)


#################
thres_list <- c(0.2,0.1,0.05,0.02,0.01,0.005,5*10^-3,10^-3, 5*10^-4, 10^-4, 5*10^-5,
                10^-5,5*10^-6, 10^-6, 5*10^-7,10^-7,5*10^-8,10^-8,0)
thres_res_1 <- find_threshold(X, Y, 1, thres_list)
idx_1 <- which(thres_res_1 == min(thres_res_1))[1]
thres_res_2 <- find_threshold(X, Y, 2, thres_list)
idx_2 <- which(thres_res_2 == min(thres_res_2))[1]
thres_res_3 <- find_threshold(X, Y, 3, thres_list)
idx_3 <- which(thres_res_3 == min(thres_res_3))[1]
thres_res_4 <- find_threshold(X, Y, 4, thres_list)
idx_4 <- which(thres_res_4 == min(thres_res_4))[1]
thres_res_5 <- find_threshold(X, Y, 5, thres_list)
idx_5 <- which(thres_res_5 == min(thres_res_5))[1]
thres_res_6 <- find_threshold(X, Y, 6, thres_list)
idx_6 <- which(thres_res_6 == min(thres_res_6))[1]
thres_res_7 <- find_threshold(X, Y, 7, thres_list)
idx_7 <- which(thres_res_7 == min(thres_res_7))[1]
thres_res_8 <- find_threshold(X, Y, 8, thres_list)
idx_8 <- which(thres_res_8 == min(thres_res_8))[1]
thres_res_9 <- find_threshold(X, Y, 9, thres_list)
idx_9 <- which(thres_res_9 == min(thres_res_9))[1]
thres_res_10 <- find_threshold(X, Y, 10, thres_list)
idx_10 <- which(thres_res_10 == min(thres_res_10))[1]
thres_res_11 <- find_threshold(X, Y, 11, thres_list)
idx_11 <- which(thres_res_11 == min(thres_res_11))[1]


# thres_val <- c(thres_res_1, thres_res_2, thres_res_3, thres_res_4, thres_res_5, thres_res_6, thres_res_7, thres_res_8, thres_res_9, thres_res_10)
# poly_order <- rep(seq(1:10), each = length(thres_list))
# thres_list_dat <- rep(thres_list, 10)
# data <- cbind(thres_list_dat, thres_val, poly_order)
# data <- data.frame(data)
# data <- subset(data, poly_order ==8)
# data <- subset(data, thres_list_dat <= thres_list[7] & thres_list_dat >= thres_list[13])
# data$poly_order <- as.factor(data$poly_order)
# p1 <- ggplot(data, aes(x=log10(thres_list_dat), y=thres_val)) 
# p1 <- p1 + geom_line()
# p1 <- p1 + ylab("Error")
# p1 <- p1 + xlab("Threshold values")
# #p1 <- p1 + scale_x_discrete(limits=-6:-1, labels =-6:-1)
# p1 <- p1 + ylim(0, 0.03)
# p1 <- p1 + theme(axis.text=element_text(size=20),
#                  axis.title=element_text(size=20,face="bold"))
# p1
thres_val <- c(thres_res_1, thres_res_2, thres_res_3, thres_res_4, thres_res_5, thres_res_6, thres_res_7, thres_res_8, thres_res_9, thres_res_10, thres_res_11)
poly_order <- rep(seq(1:11), each = length(thres_list))
thres_list_dat <- rep(thres_list, 11)
data <- cbind(thres_list_dat, thres_val, poly_order)
data <- data.frame(data)
data_1 <- subset(data, poly_order ==3 | poly_order ==4)
#data <- subset(data, thres_list_dat >=0 & thres_list_dat < 0.02)
data_1$poly_order <- as.factor(data_1$poly_order)

p1 <- ggplot(data_1, aes(x=log10(thres_list_dat),y=thres_val, group = poly_order)) 
p1 <- p1 + geom_line(aes(color = poly_order, linetype = poly_order))
p1 <- p1 + ylab("Cross-validation Error")
p1 <- p1 + xlab("Threshold values (log)")
p1 <- p1 + ylim(4, 11)
p1 <- p1 + theme(axis.text=element_text(size=20),
                 axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=20),
                 legend.justification=c(0,1),legend.position=c(0,1),legend.title=element_blank())
p1 <- p1 +  scale_linetype_manual(values = c(3, 2,1 ),
                                  labels = c("p = 3", "p = 4")) 
p1 <- p1 + geom_point(x= log10(thres_list[idx_3]), y = 8.167529, color='blue')
p1 <- p1 + geom_point(x= log10(thres_list[idx_4]), y = 4.737542, color='red')
#p1 <- p1 +  geom_vline(xintercept=log10(thres_list[idx_4]), linetype="dashed", 
#                       color = "red", size=1)
#p1 <- p1 +  geom_vline(xintercept=log10(thres_list[idx_3]), linetype="dashed", 
#                       color = "blue", size=1)
p1 <- p1 + scale_color_manual(values = c("blue", "brown", "black"),
                              labels = c("p = 3", "p = 4")) 
p1


thres_val <- c(thres_res_1, thres_res_2, thres_res_3, thres_res_4, thres_res_5, thres_res_6, thres_res_7, thres_res_8, thres_res_9, thres_res_10, thres_res_11)
poly_order <- rep(seq(1:11), each = length(thres_list))
thres_list_dat <- rep(thres_list, 11)
data <- cbind(thres_list_dat, thres_val, poly_order)
data <- data.frame(data)
data_2 <- subset(data, poly_order ==8 | poly_order ==9)
data_2 <- subset(data_2, thres_list_dat > 2*10^-5 & thres_list_dat < 0.03)
data_2$poly_order <- as.factor(data_2$poly_order)

p2 <- ggplot(data_2, aes(x=log10(thres_list_dat),y=thres_val, group = poly_order)) 
p2 <- p2 + geom_line(aes(color = poly_order, linetype = poly_order))
p2 <- p2 + ylab("Cross-validation Error")
p2 <- p2 + xlab("Threshold values (log)")
p2 <- p2 + ylim(0, 0.3)
p2 <- p2 + theme(axis.text=element_text(size=20),
                 axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=20),
                 legend.justification=c(0,1),legend.position=c(0,1),legend.title=element_blank())
p2 <- p2 +  scale_linetype_manual(values = c(3, 2,1 ),
                                  labels = c("p = 8", "p = 9")) 
p2 <- p2 + geom_point(x= log10(0.001), y = 0.006400592 , color='blue',size = 3)
p2 <- p2 + geom_point(x= log10(0.001), y = 0.004723677 , color='red')

#p2 <- p2 +  geom_vline(xintercept=log10(0.001), linetype="dashed", 
#                       color = "red", size=1)
p2 <- p2 + scale_color_manual(values = c("blue", "brown", "black"),
                              labels = c("p = 8", "p = 9")) 
p2


# thres_val <- c(thres_res_1, thres_res_2, thres_res_3, thres_res_4, thres_res_5, thres_res_6, thres_res_7, thres_res_8, thres_res_9, thres_res_10, thres_res_11)
# poly_order <- rep(seq(1:11), each = length(thres_list))
# thres_list_dat <- rep(thres_list, 11)
# data <- cbind(thres_list_dat, thres_val, poly_order)
# data <- data.frame(data)
# data <- subset(data, poly_order >=10)
# data <- subset(data, thres_list_dat > 10^-8 & thres_list_dat < 0.05)
# data$poly_order <- as.factor(data$poly_order)
# 
# p1 <- ggplot(data, aes(x=log10(thres_list_dat),y=thres_val, group = poly_order)) 
# p1 <- p1 + geom_line(aes(color = poly_order, linetype = poly_order))
# p1 <- p1 + ylab("Error")
# p1 <- p1 + xlab("Threshold values (log)")
# #p1 <- p1 + ylim(0, 0.03)
# p1 <- p1 + theme(axis.text=element_text(size=20),
#                  axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=20),
#                  legend.justification=c(0,1),legend.position=c(0,1),legend.title=element_blank())
# p1 <- p1 +  scale_linetype_manual(values = c(3, 2,1 ),
#                                    labels = c("p = 8", "p = 9", "p = 10")) 
# p1 <- p1 + scale_color_manual(values = c("blue", "brown", "black"),
#                               labels = c("p = 8", "p = 9", "p = 10")) 
# p1

############# Plot standard deviation error graph 
var = 49/8 + 0.1*pi^4/5 + 0.01*pi^8/18 + 0.5
sd_true <- sqrt(var)
sd_true
LAR_mom_approx <- function(X, Y, p){
  X_feature <- create_X_feature(X, p)
  X_feature <- scale(X_feature)
  X_ortho <- gramschmidt(X_feature)$Q
  #print(rankMatrix(X_ortho))
  lar.1 <- lars(y=Y, x= X_ortho, type="lar", normalize = FALSE, trace=FALSE, intercept = TRUE, Gram = FALSE, use.Gram = FALSE, max.steps = 50)
  cv.lar.1 <- cv.lars(y=Y, x= X_ortho, type="lar", K=5, mode = "step", intercept = TRUE, Gram = FALSE, use.Gram = FALSE, max.steps = 50)
  s_cv <- which.min(cv.lar.1$cv) 
  coef_lar <- coef(lar.1, s=s_cv, mode="step") 
  index <- coef_lar !=0 
  m1 <- lm(Y~1 + X_ortho[,index])
  mean_approx_lar <- m1$coefficients[1]
  sd_approx_lar <- sqrt(sum(m1$coefficients[-1]^2))
  return(list('mean'=mean_approx_lar, 'sd'=sd_approx_lar))
}




# LAR_mom_approx <- function(X, Y, p){
#   X_feature <- create_X_feature(X, p)
#   X_feature <- scale(X_feature)
#   X_ortho <- gramschmidt(X_feature)$Q
#   lar.1 <- lars(y=Y, x= X_ortho, type="lasso", trace=TRUE)
#   #cv.lar.1 <- cv.lars(y=Y, x= X_ortho, type="lasso", K=5, mode = "step")
#   a <- summary(lar.1)
#   s_cv <- which.min(a$Cp)
#   #s_cv <- which.min(cv.lar.1$cv)
#   coef_lar <- coef(lar.1, s=85, mode="step")
#   index <- coef_lar !=0
#   m1 <- lm(Y~1 + X_ortho[,index])
#   mean_approx_lar <- m1$coefficients[1]
#   sd_approx_lar <- sqrt(sum(m1$coefficients[-1]^2))
#   return(list('mean'=mean_approx_lar, 'sd'=sd_approx_lar))
# }


data <- data_create(3, 200, LHS = TRUE)
X <- data$X
Y <- data$Y

data_lar <- data_create(3, 400, LHS = TRUE)
X_lar <- data_lar$X
Y_lar <- data_lar$Y

l2_error_p3 <- rep(0, 11)
lar_error_p3_small <- rep(0, 11)
lar_error_p3_large <- rep(0, 11)


for(i in 1:11){
  print(i)
  l2_temp <- rep(0,50)
  lar_temp_small <- rep(0, 50)
  lar_temp_large <- rep(0, 50)
  for(j in 1:50){
    data <- data_create(j, 200, LHS = TRUE)
    X <- data$X
    Y <- data$Y
    data_lar <- data_create(j, 400, LHS = TRUE)
    X_lar <- data_lar$X
    Y_lar <- data_lar$Y
    thres_res <- find_threshold(X, Y, i, thres_list)
    idx <- which(thres_res == min(thres_res))[1]
    sd_approx<- mom_approxi_propose(X, Y, i, thres_list[idx])$sd_approx
    l2_temp[j] <- abs(sd_approx - sd_true)/sd_true
    
    lar_approx_small <- LAR_mom_approx(X, Y, i)$sd
    lar_temp_small[j] <- abs(lar_approx_small - sd_true)/sd_true 
    
    lar_approx_large <- LAR_mom_approx(X_lar, Y_lar, i)$sd
    lar_temp_large[j] <- abs(lar_approx_large - sd_true)/sd_true 
    
    
  }
  l2_error_p3[i] <- mean(l2_temp)
  lar_error_p3_small[i] <- mean(lar_temp_small)
  lar_error_p3_large[i] <- mean(lar_temp_large)
  #print(sd_approx)
}

error <- c(l2_error_p3[1:10], lar_error_p3_small[1:10])
data_error <- cbind(rep(seq(1:10), times=2), error, rep(c("proposed method", "LAR-based"), each = 10))
data_error <- data.frame(data_error)
colnames(data_error) <- c('poly_order', 'error', 'type')
data_error$poly_order<- as.numeric(as.character(data_error$poly_order))
data_error$error <- as.numeric(as.character(data_error$error))

p3 <- ggplot(data_error, aes(x=poly_order, y= error, group = type))
p3 <- p3 + geom_line(aes(color = type, linetype = type))
p3 <- p3 + ylab("Error")
p3 <- p3 + xlab("Polynomial order")
p3 <- p3 + scale_x_discrete(limits=1:10, labels =1:10)
p3 <- p3 + ylim(0, 0.6)
p3 <- p3 + theme(axis.text=element_text(size=20),
                 axis.title=element_text(size=20,face="bold"),
                 legend.text=element_text(size=20),
                 legend.justification=c(1,1),legend.position=c(1,1),legend.title=element_blank())
p3




error <- c(l2_error_p[1:10], lar_error_p[1:10])
data_error <- cbind(rep(seq(1:10), times=2), error, rep(c("proposed method", "LAR-based"), each = 10))
data_error <- data.frame(data_error)
colnames(data_error) <- c('poly_order', 'error', 'type')
data_error$poly_order<- as.numeric(as.character(data_error$poly_order))
data_error$error <- as.numeric(as.character(data_error$error))

p3 <- ggplot(data_error, aes(x=poly_order, y= error, group = type))
p3 <- p3 + geom_line(aes(color = type, linetype = type))
p3 <- p3 + ylab("Error")
p3 <- p3 + xlab("Polynomial order")
p3 <- p3 + scale_x_discrete(limits=1:10, labels =1:10)
p3 <- p3 + ylim(0, 0.6)
p3 <- p3 + theme(axis.text=element_text(size=20),
                 axis.title=element_text(size=20,face="bold"),
                 legend.text=element_text(size=20),
                 legend.justification=c(1,1),legend.position=c(1,1),legend.title=element_blank())
p3




n_list <- c(100, 200, 300, 400, 500, 600, 700, 800)
thres_list <- c(0.2,0.1,0.05,0.02,0.01,0.005,5*10^-3,10^-3, 5*10^-4, 10^-4, 5*10^-5,
                10^-5,5*10^-6, 10^-6, 5*10^-7,10^-7,5*10^-8,10^-8,0)

l2_error_sample <- rep(0, length(n_list))
lar_error_sample <- rep(0, length(n_list))
l2_std <- rep(0, length(n_list))
lar_std <- rep(0, length(n_list))
for(i in 1:length(n_list)){
  print(i)
  l2_temp <- rep(0, 50)
  lar_temp <- rep(0, 50)
  n <- n_list[i]
  for(j in 1:50){
    #print(j)
    data <- data_create(j, n, LHS = TRUE)
    X <- data$X
    Y <- data$Y
    # if(n <400){
    #   thres_res_11 <- find_threshold(X, Y, 8, thres_list)
    #   idx_11 <- which(thres_res_11 == min(thres_res_11))[1]
    #   sd_approx<- mom_approxi_propose(X, Y, 8, thres_list[idx])$sd_approx
    #   l2_temp[j] <- abs(sd_approx - sd_true)/sd_true
    # }
    # else{
    thres_res_11 <- find_threshold(X, Y, 8, thres_list)
    idx_11 <- which(thres_res_11 == min(thres_res_11))[1]
    sd_approx<- mom_approxi_propose(X, Y, 8, thres_list[idx])$sd_approx
    l2_temp[j] <- abs(sd_approx - sd_true)/sd_true
    lar_approx <- LAR_mom_approx(X, Y, i)$sd
    lar_temp[j] <- abs(lar_approx - sd_true)/sd_true
    #}
    l2_error_sample[i] <- mean(l2_temp)
    lar_error_sample[i] <- mean(lar_temp)
    l2_std[i] <- sd(l2_temp)
    lar_std[i] <- sd(lar_temp)
  }
}

#lar_error_sample[1] <- NA
#lar_error_sample[2] <- NA
#lar_error_sample[3] <- NA

error <- c(l2_error_sample, lar_error_sample)
data_error <- cbind(rep(n_list, times=2), error, rep(c("proposed method", "LAR-based"), each = 8))
data_error <- data.frame(data_error)
colnames(data_error) <- c('sample_size', 'error', 'type')
data_error$sample_size<- as.numeric(as.character(data_error$sample_size))
data_error$error <- as.numeric(as.character(data_error$error))

p4 <- ggplot(data_error, aes(x=sample_size, y= error, group = type))
p4 <- p4 + geom_line(aes(color = type, linetype = type))
p4 <- p4 + ylab("Error")
p4 <- p4 + xlab("Sample size")
p4 <- p4 + scale_x_discrete(limits=n_list, labels =n_list)
p4 <- p4 + ylim(0, 0.2)
p4 <- p4 + theme(axis.text=element_text(size=20),
                 axis.title=element_text(size=20,face="bold"),
                 legend.text=element_text(size=20),
                 legend.justification=c(1,1),legend.position=c(1,1),legend.title=element_blank())
p4





data_lar <- data_create(3, 1000, LHS = TRUE)
X_lar <- data_lar$X
Y_lar <- data_lar$Y

l2_time <- rep(0, 13)
lar_time <- rep(0, 13)
l2_time_std <- rep(0, 13)
lar_time_std <- rep(0, 13)
for(i in 1:13){
  print(i)
  l2_temp <- rep(0, 50)
  lar_temp <- rep(0, 50)
  for(j in 1:50){
    thres_res <- find_threshold(X, Y, i, thres_list)
    idx <- which(thres_res == min(thres_res))[1]
    l2_start <- Sys.time()
    X_feature <- create_X_feature(X, i)
    res <- sparse_GS(X_feature, Y, thres_list[idx])
    l2_end <- Sys.time()
    
    lar_start <- Sys.time()
    X_feature <- create_X_feature(X, i)
    X_feature <- scale(X_feature)
    X_ortho <- gramschmidt(X_feature)$Q
    lar_end <- Sys.time()
    l2_temp[j] <- l2_end - l2_start
    lar_temp[j] <- lar_end - lar_start 
  }

  #print(sd_approx)
  l2_time[i] <- mean(l2_temp)
  lar_time[i] <- mean(lar_temp)
  l2_time_std[i] <- sd(l2_temp)
  lar_time_std[i] <- sd(lar_temp)
}



time <- c(l2_time, lar_time)
data_time <- cbind(rep(seq(1:13), times=2), time, rep(c("proposed method", "LAR-based"), each = 13))
data_time <- data.frame(data_time)
colnames(data_time) <- c('polynomial_order', 'time', 'type')
data_time$polynomial_order<- as.numeric(as.character(data_time$polynomial_order))
data_time$time <- as.numeric(as.character(data_time$time))

p5 <- ggplot(data_time, aes(x=polynomial_order, y= time, group = type))
p5 <- p5 + geom_line(aes(color = type, linetype = type))
p5 <- p5 + ylab("Time")
p5 <- p5 + xlab("Polynomial order")
p5 <- p5 + scale_x_discrete(limits=seq(1:13), labels =seq(1:13))
p5 <- p5 + ylim(0, 2)
p5 <- p5 + theme(axis.text=element_text(size=20),
                 axis.title=element_text(size=20,face="bold"),
                 legend.text=element_text(size=20),
                 legend.justification=c(0,1),legend.position=c(0,1),legend.title=element_blank())
p5










thres_val <- c(thres_res_1, thres_res_2, thres_res_3, thres_res_4, thres_res_5, thres_res_6, thres_res_7, thres_res_8, thres_res_9, thres_res_10)
poly_order <- rep(seq(1:10), each = length(thres_list))
thres_list_dat <- rep(thres_list, 10)
data <- cbind(thres_list_dat, thres_val, poly_order)
data <- data.frame(data)
data <- subset(data, poly_order ==8)
data <- subset(data, thres_list_dat <= thres_list[7] & thres_list_dat >= thres_list[13])
data$poly_order <- as.factor(data$poly_order)
ggplot(data, aes(x=log10(thres_list_dat), y=thres_val, group = poly_order)) + geom_line(aes(color=poly_order, linetype =poly_order))


#thres_res_9
thres_res_10 <- find_threshold(X, Y, 10, thres_list)
idx_10 <- which(thres_res_10 == min(thres_res_10))[1]

# thres_res_12 <- find_threshold(X, Y, 12, thres_list)

sd_approx <- matrix(0, nrow = 100, ncol = 10)
mean_approx <- matrix(0, nrow = 100, ncol = 10)
for(i in 1:100){
  data <- data_create(i, 1000, LHS = TRUE)
  X <- data$X
  Y <- data$Y
  sd_approx[i,1] <- mom_approxi_propose(X, Y, 1, thres_list[idx_1])$sd_approx
  sd_approx[i,2] <- mom_approxi_propose(X, Y, 2, thres_list[idx_2])$sd_approx
  sd_approx[i,3] <- mom_approxi_propose(X, Y, 3, thres_list[idx_3])$sd_approx
  sd_approx[i,4] <- mom_approxi_propose(X, Y, 4, thres_list[idx_4])$sd_approx
  sd_approx[i,5] <- mom_approxi_propose(X, Y, 5, thres_list[idx_5])$sd_approx
  sd_approx[i,6] <- mom_approxi_propose(X, Y, 6, thres_list[idx_6])$sd_approx
  sd_approx[i,7] <- mom_approxi_propose(X, Y, 7, thres_list[idx_7])$sd_approx
  sd_approx[i,8] <- mom_approxi_propose(X, Y, 8, thres_list[idx_8])$sd_approx
  sd_approx[i,9] <- mom_approxi_propose(X, Y, 9, thres_list[idx_9])$sd_approx
  sd_approx[i,10] <- mom_approxi_propose(X, Y, 10, thres_list[idx_10])$sd_approx
}

apply(sd_approx, 2, mean)
apply(sd_approx, 2, sd)

mom_approxi_propose(X, Y, 1, thres_list[idx_1])
mom_approxi_propose(X, Y, 2, thres_list[idx_2])
mom_approxi_propose(X, Y, 3, thres_list[idx_3])
mom_approxi_propose(X, Y, 4, thres_list[idx_4])
mom_approxi_propose(X, Y, 5, thres_list[idx_5])
mom_approxi_propose(X, Y, 6, thres_list[idx_6])
mom_approxi_propose(X, Y, 7, thres_list[idx_7])
mom_approxi_propose(X, Y, 8, thres_list[idx_8])
mom_approxi_propose(X, Y, 9, thres_list[idx_9])




var = 49/8 + 0.1*pi^4/5 + 0.01*pi^8/18 + 0.5
sd_true <- sqrt(var)

####LARs based moment approximation 
LAR_mom_approx <- function(X, Y, p){
  X_feature <- create_X_feature(X, p)
  X_feature <- scale(X_feature)
  X_ortho <- gramschmidt(X_feature)$Q
  lar.1 <- lars(y=Y, x= X_ortho, type="lasso", trace=TRUE)
  #a <- summary(lar.1)
  #coef_lar <- coef(lar.1, s=which.min(a$Cp), mode="step") 
  
  # Print out coefficients at optimal s.
  cv.lar.1 <- cv.lars(y=Y, x= X_ortho, type="lasso", K=5, mode = "step")
  s_cv <- which.min(cv.lar.1$cv)
  coef_lar <- coef(lar.1, s=s_cv, mode="step") 
  index <- coef_lar !=0 
  m1 <- lm(Y~1 + X_ortho[,index])
  mean_approx_lar <- m1$coefficients[1]
  sd_approx_lar <- sqrt(sum(m1$coefficients[-1]^2))
  return(list('mean'=mean_approx_lar, 'sd'=sd_approx_lar))
}

LAR_mom_approx(X, Y, 6)



### analytic mean 3.5  analytic std 3.7208 
data <- data_create(3, 1000, LHS = TRUE)
X <- data$X
Y <- data$Y
res <- mom_approxi_propose(X, Y, 10, 0)
res$mean_approx
res$sd_approx




### plot mean and variance convergency plot for both algorithms refer to aPCE paper in RESS

## https://arxiv.org/pdf/1907.03933.pdf





################### LAR related codes for reference codes 
data <- data_create(3, 1000, LHS = TRUE)
X <- data$X
Y <- data$Y
X_feature <- create_X_feature(X, 10)
X_ortho <- gramschmidt(X_feature)$Q
res <- cv.lars(X_ortho, Y, K = 5, type = 'forward.stagewise', mode ='step')
res$cv.error
res$mode
res$cv
coef(res)


cvlas <- cv.lars(X_ortho, Y, type="lasso")
cvlas$cv


object <- lars(X_ortho, Y,type="lasso")
coef(object)


las.coef <- predict.lars(las, type="coefficients", mode="fraction", s=frac)

cvlas <- cv.lars(X_ortho, Y, type="lasso")
cvlas


las <- lars(X_ortho, Y, type="lasso")
las
plot(las, plottype="coefficients")


frac <- cvlas$fraction[which.min(cvlas$cv)]
frac
las.coef <- predict.lars(las, type="coefficients", mode="fraction", s=frac)
las.coef



# First half of data 
data <- data_create(3, 1000, LHS = TRUE)
X <- data$X
Y <- data$Y
X_feature <- create_X_feature(X, 6)
X_ortho <- gramschmidt(X_feature)$Q
lar.1 <- lars(y=Y, x= X_ortho, type="lasso", trace=TRUE)
summary(lar.1)
#plot(lar.1) # Plots coefficient path
#coef(lar.1) # Lists out coefficients for each step in path
# cv.lars() uses crossvalidation to estimate optimal position in path
cv.lar.1 <- cv.lars(y=Y, x= X_ortho, type="lasso", K=5)
cv.lar.1

# Use the best Cp value to find best model: 
a <- summary(lar.1)
# Print out coefficients at optimal s.
s_cv <- which.min(cv.lar.1$cv)
coef_lar1 <- coef(lar.1, s=which.min(a$Cp), mode="step") 
sqrt(sum(coef_lar1[-1]^2))

coef_lar <- coef(lar.1, s=s_cv, mode="step") 
coef_lar[1]
sqrt(sum(coef_lar^2))
mean(Y)
sd(Y)
