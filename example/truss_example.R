library(VineCopula)
library(stats)
library(evd)
library(CDVine)
library(VineCopula)
library(stats)
library(evd)
library(CDVine)
data_create_3 <- function(seed, size){
  set.seed(seed)
  X1 <- create_lognorm(size, 2.1*10^11, 2.1*10^10)
  X2 <- create_lognorm(size, 2.1*10^11, 2.1*10^10)
  X3 <- create_lognorm(size, 2.0*10^-3, 2.0*10^-4)
  X4 <- create_lognorm(size, 1.0*10^-3, 1.0*10^-4)
  d <- 6 
  dd <- d*(d-1)/2
  fam <- rep(4, dd)
  par <- c(1.1,1.1,1.1,1.1,1.1, 
           1,1,1,1,1,
           1,1,1,1,1)
  U = CDVineSim(size,fam,par,type=1)
  mu <- 5*10^4
  sigma <- 0.15*mu
  beta <- sqrt(6)*sigma/pi
  alpha <- mu - 0.5772*beta
  X5 <- create_gumbel(alpha, beta, U[,1])
  X6 <- create_gumbel(alpha, beta, U[,2])  
  X7 <- create_gumbel(alpha, beta, U[,3]) 
  X8 <- create_gumbel(alpha, beta, U[,4])
  X9 <- create_gumbel(alpha, beta, U[,5])
  X10 <- create_gumbel(alpha, beta, U[,6])
  
  X <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
  
  X1 <- (X1 - 2.1*10^11)/(2.1*10^10)
  X2 <- (X2 - 2.1*10^11)/(2.1*10^10)
  X3 <- (X3 - 2.0*10^-3)/(2.0*10^-4)
  X4 <- (X4 - 1.0*10^-3)/(1.0*10^-4)
  X5 <- (X5 - 5*10^4)/(0.15*mu)
  X6 <- (X6 - 5*10^4)/(0.15*mu)
  X7 <- (X7 - 5*10^4)/(0.15*mu)
  X8 <- (X8 - 5*10^4)/(0.15*mu)
  X9 <- (X9 - 5*10^4)/(0.15*mu)
  X10 <- (X10 - 5*10^4)/(0.15*mu)
  
  
  Y <- 2.8070 + 1.2598*X1 + 0.2147*X2 + 1.2559*X3 + 0.2133*X4 - 0.1510*X5 -0.4238*X6-0.6100*X7-0.6100*X8 - 0.4238*X9 - 0.1510*X10 - 0.1978*X1^2 - 
    0.0362*X2^2 -0.2016*X3^2 - 0.0346*X4^2+0.0023*X5^2 + 0.0008*X6^2 + 0.0036*X7^2 + 0.0036*X8^2 + 0.0008*X9^2 + 0.0023*X10^2 - 0.0042*X1*X2 -
    0.3022*X1*X3 - 0.0110*X1*X4 + 0.0381*X1*X5 + 0.0871*X1*X6 + 0.1232*X1*X7 + 0.1232*X1*X8 + 
    0.0871*X1*X9 + 0.0346*X1*X10 + 0.0041*X2*X3 + 0.0110*X3*X4 + 0.0261*X3*X5 + 0.0831*X3*X6+ 
    0.1172*X3*X7 + 0.1172*X3*X8 + 0.0832*X3*X9 + 0.0296*X3*X10
  #X <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
  return(list('X'=X, 'Y'=Y))
}




create_gumbel<- function(loc, scale, Z){
  res <- rep(0, length(Z))
  for(i in 1:length(Z)){
    res[i] <- qgumbel(Z[i], loc, scale)
  }
  return(res)
}


create_lognorm <- function(size, m, s){
  location <- log(m^2 / sqrt(s^2 + m^2))
  shape <- sqrt(log(1 + (s^2 / m^2)))
  res <- rlnorm(size, location, shape)
  return(res)
}



sd_true_list <- rep(0, 100)
for(j in 1:100){
  data_temp <- data_create_3(j,10^5)  
  sd_true_list[j] <- sd(data_temp$Y)
}

mean(sd_true_list)
sd(sd_true_list)


res <- function(seed, size, thres_list, sd_true){
  data <- data_create_3(seed, size)
  X <- data$X
  Y <- data$Y
  thres_res <- find_threshold(X, Y, 2, thres_list)
  idx <- which(thres_res == min(thres_res))[1]
  print(thres_list[idx])
  X_feature <- create_X_feature(X, 2)
  res <- sparse_GS(X_feature, Y, thres_list[idx])
  num_l2 <- max(res$select_idx)
  m1 <- lm(Y ~res$ortho_vec[,-1])
  coef <- m1$coefficients[!is.na(m1$coefficients)]
  sd_approx <- sqrt(sum(coef[-1]^2))
  l2_error <- abs(sd_approx - sd_true)/sd_true
  
  X_feature <- create_X_feature(X, 2)
  X_feature <- scale(X_feature)
  X_ortho <- gramschmidt(X_feature)$Q
  lar.1 <- lars(y=Y, x= X_ortho, type="lar", normalize = FALSE, trace=FALSE, intercept = TRUE, Gram = FALSE, use.Gram = FALSE, max.steps = size_1)
  cv.lar.1 <- cv.lars(y=Y, x= X_ortho, type="lar", K=5, mode = "step", intercept = TRUE, Gram = FALSE, use.Gram = FALSE, max.steps = size_1)
  s_cv <- which.min(cv.lar.1$cv) 
  coef_lar <- coef(lar.1, s=s_cv, mode="step") 
  index <- coef_lar !=0 
  if(s_cv == 1){
    m1 <- lm(Y~1) 
  }
  else{
    m1 <- lm(Y~1 + X_ortho[,index])
  }
  num_lar <- sum(index) -1
  sd_approx_lar <- sqrt(sum(m1$coefficients[-1]^2))
  lar_error <- abs(sd_approx_lar - sd_true)/sd_true
  
  return(list('num_l2' = num_l2, 'num_lar' = num_lar, 'l2_error' = l2_error, 'lar_error'=lar_error,
              'l2_approx' =sd_approx, 'lar_approx' =sd_approx_lar))
}

data <- data_create_3(3, 10^6)
cor(data$X[,5:10])


thres_list <- c(10^-2,10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8, 0)


l2_error_small_samp <- rep(0, 50)
lar_error_small_samp <- rep(0, 50)
l2_approx_small <- rep(0, 50)
lar_approx_small <- rep(0,50)
for(i in 1:50){
  print(i)
  temp_res <- res(i, 50, thres_list, sd_true)
  l2_error_small_samp[i] <- temp_res$l2_error
  lar_error_small_samp[i] <- temp_res$lar_error
  l2_approx_small[i] <- temp_res$l2_approx
  lar_approx_small[i] <- temp_res$lar_approx
}
mean(l2_approx_small)
mean(lar_approx_small)

sd(l2_approx_small)/sqrt(50)
sd(lar_approx_small)/sqrt(50)





l2_error_large_samp <- rep(0, 50)
lar_error_large_samp <- rep(0, 50)
l2_approx_large <- rep(0, 50)
lar_approx_large <- rep(0,50)
for(i in 1:50){
  print(i)
  temp_res <- res(i, 100,  thres_list, sd_true)
  l2_error_large_samp[i] <- temp_res$l2_error
  lar_error_large_samp[i] <- temp_res$lar_error
  l2_approx_large[i] <- temp_res$l2_approx
  lar_approx_large[i] <- temp_res$lar_approx
}
mean(l2_approx_large)
mean(lar_approx_large)

sd(l2_approx_large)/sqrt(50)
sd(lar_approx_large)/sqrt(50)


data <- data_create_3(3, 1000)
X <- data$X
Y <- data$Y
thres_list <- c(10^-12,3*10^-13, 10^-14, 10^-15,0)
thres_res_3 <- find_threshold(X, Y, 3, thres_list)
thres_res_3


res <- mom_approxi_propose(X, Y, 3, 10^-15)
res

sd_true <- sd(data_create_3(3, 10^5)$Y)*sqrt((10^5-1)/10^5) 

sample_size <- c(100, 500, 1000, 1500, 2000, 2500, 3000)
l2_error <- rep(0, length(sample_size))
lar_error <- rep(0, length(sample_size))
for(i in 1:length(sample_size)){
  print(3)
  data <- data_create_3(i, sample_size[i])
  X <- data$X
  Y <- data$Y
  res <- mom_approxi_propose(X, Y, 3, 10^-15)
  l2_error[i] <- abs(res$sd_approx- sd_true)/sd_true
  lar_error[i] <- abs(LAR_mom_approx(X,Y,3)$sd -sd_true)/sd_true
}


error <- c(l2_error, lar_error)
data_error <- cbind(rep(sample_size, times=2), error, rep(c("proposed method", "LAR-based"), each = length(sample_size)))
data_error <- data.frame(data_error)
colnames(data_error) <- c('sample_size', 'error', 'type')
data_error$sample_size<- as.numeric(as.character(data_error$sample_size))
data_error$error <- as.numeric(as.character(data_error$error))

p2 <- ggplot(data_error, aes(x=sample_size, y= error, group = type))
p2 <- p2 + geom_line(aes(color = type, linetype = type))
p2 <- p2 + ylab("Error")
p2 <- p2 + xlab("Sample size")
#p2 <- p2 + scale_x_discrete(limits=1:length(sample_size), labels =sample_size)
p2 <- p2 + ylim(0, 0.15)
p2 <- p2 + theme(axis.text=element_text(size=20),
                 axis.title=element_text(size=20,face="bold"),
                 legend.text=element_text(size=20),
                 legend.justification=c(1,1),legend.position=c(1,1),legend.title=element_blank())
p2





### plot mean and variance convergency plot for both algorithms 














obj <- BiCop(family = 2, par = 0.4, par2 = 6)
obj
## see the object's content or a summary
str(obj)
summary(obj)
## a selection of functions that can be used with BiCop objects
simdata <- BiCopSim(300, obj)

mu <- 5*10**4
sigma <- 0.15*mu
beta <- sqrt(6)*sigma/pi
alpha <- mu - 0.5772*beta


obj <- BiCop(family = 4, par= 1.1)
simdata <- BiCopSim(300, obj)

d <- 6 
dd <- d*(d-1)/2
fam <- rep(4, dd)
par <- c(1.1,1.1,1.1,1.1,1.1, 
           1,1,1,1,1,
           1,1,1,1,1)
N <- 100
U1 = CDVineSim(N,fam,par,type=1)
head(U1)


d = 6
dd = d*(d-1)/2
fam1 = rep(1,dd)
par1 = c(0.2,0.69,0.73,0.22,-0.09,0.51,0.32,0.01,0.82,0.01,
         -0.2,-0.32,-0.19,-0.17,-0.06)
N = 100
U1 = CDVineSim(N,fam1,par1,type=1)
head(U1)

