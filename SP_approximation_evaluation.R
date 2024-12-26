library(tidyverse)
library(cubature)
library(akima)
library(nor1mix)
library(LaplacesDemon)

library(foreach)
library(doParallel)
library(doSNOW)

load("C:/Zehang_workspace/Approximation_Supply_Demand/day_ahead_original_supply_curves.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/day_ahead_original_demand_curves.RData")


####################################################

#---- 0. Defining some useful functions----                    

####################################################



##### 0.1 distance function ####
stepwise_distance_w_plot <- function(p1, q1, #prices and quantities for curve 1
                                     p2, q2, #prices and quantities for curve 2
                                     l, #l of the distance
                                     w_type = c("none", "max_price", "user_wfun"), #type of w(p)
                                     maxprice, # only for max_price. Maximum price considered.
                                     user_fun, # only for user_fun. Function that integrates up to one.
                                     plot = TRUE){
  
  # Define the two supply functions as R functions
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  # Define the weighting function
  # Option 1: W(p) = 1. Last jumps are not considered
  if(w_type == "none"){
    
    # Define a grid of evaluation points
    p <- unique(sort(c(0, p1, p2))) #, Inf)))
    grid <- (p[-1]+p[-length(p)])/2
    
    # Define the weighting function W(p) = 1
    swfun  <- stepfun(-1, c(1, 1), f = 0)
    
    w <- sapply(1:c(length(p)-1), function(x) integrate(swfun, 
                                                        upper = p[-1][x], 
                                                        lower = p[-length(p)][x])$value)
    
    step_dist <- sum(abs(sfun1(grid)-sfun2(grid))^l*w)
  }
  # Option 2: W(p) = 1 until a p < maxprice
  if(w_type == "max_price"){
    if(missing(maxprice)){maxprice <- 10000}
    
    # Define a grid of evaluation points
    p <- unique(sort(c(0, p1, p2, maxprice)))
    grid <- (p[-1]+p[-length(p)])/2
    
    swfun  <- stepfun(maxprice, c(1, 0), f = 0)
    
    w <- sapply(1:c(length(p)-1), function(x) integrate(swfun, 
                                                        upper = p[-1][x], 
                                                        lower = p[-length(p)][x])$value)
    
    step_dist <- sum(abs(sfun1(grid)-sfun2(grid))^l*w)
  }
  
  if(w_type == "user_wfun"){
    swfun <- user_fun
    # Define a grid of evaluation points
    p <- unique(sort(c(0, p1, p2, Inf))) # W(t) goes to 0
    grid <- (p[-1]+p[-length(p)])/2
    
    w <- sapply(1:c(length(p)-1), function(x) integrate(swfun, 
                                                        upper = p[-1][x], 
                                                        lower = p[-length(p)][x])$value)
    
    step_dist <- sum(abs(sfun1(grid)-sfun2(grid))^l*w)
  }
  
  if(plot == TRUE & w_type %in% c("none", "max_price", "user_wfun")){
    layout(matrix(c(1, 1, 1, 1, 2, 2), nrow = 3, ncol = 2, byrow = T))
    
    if(missing(maxprice)){maxprice = 0}
    plot(sfun1, xlim = c(0, max(p1, p2, maxprice)), ylim = c(0, max(q1, q2)),
         ylab = "Quantity", xlab = "Price", 
         main = paste0("l", l," Stepwise Distance = ", step_dist))
    plot(sfun2, col = "blue",  xlim = c(0, max(p1, p2, maxprice)), add = TRUE)
    points(c(0, p1), q1, pch = 16)
    points(c(0, p2), q2, pch = 16, col = "blue")
    
    plot(swfun, xlim = c(0, max(p1, p2, maxprice)),
         ylab = "Weight", xlab = "Price", 
         main = paste0("Weigthing function: ", w_type))
  }
  
  return(step_dist)
}


##### 0.2 weight function being used in the SUPPLY CURVE approximation procedure ####

load("C:/Zehang_workspace/Approximation_Supply_Demand/norMixmodel_supply_p_scaled.RData")
supply.mixture.gaussian = function(x) dnorMix(norMixmodel_supply_p_scaled,x)
# supply_p_scaled = day_ahead_supply$price/max(day_ahead_supply$price)
# norMixmodel_supply_p_scaled <-  norMixMLE(supply_p_scaled,m=4)
# plot(density(supply_p_scaled),col=3)
# lines(norMixmodel_supply_p_scaled)


##### 0.3 weight function being used in the DEMAND CURVE approximation procedure ####
load("C:/Zehang_workspace/Approximation_Supply_Demand/norMixmodel_demand_p_scaled.RData")
demand.mixture.gaussian = function(x) dnorMix(norMixmodel_demand_p_scaled,x)
# demand_p_scaled = day_ahead_demand$price / max(day_ahead_demand$price)
# norMixmodel_demand_p_scaled = norMixMLE(demand_p_scaled,m=4)
# plot(density(demand_p_scaled),col=3)
# lines(norMixmodel_demand_p_scaled)


# ##### 0.4 nodes selection function ####
# # margina distribution
# f_EmprDis= function(p,dat) mean(p<=dat)
# 
# # conditional distribution
# f_EmprBivDis_Q = function(p,Q,dat) {
#   dat = dat[dat[,2]>=Q,]
#   return(f_EmprDis(p,dat[,1]))
# } 


##### 0.5 function for dyadic nodes selecting
fapprox =  function(x,y,eval) {
  x = x[-1]
  sf = stepfun(x, y, f = 0) 
  return(sf(eval))}

f_interval_approx_error = function(cp,cq, # the p and q of all jumps
                                   l_boundary, # the left boundary of interval
                                   r_boundary, # the right boundary of interval
                                   f_weight,
                                   l ){# L1 or L2 norm
  
  I_cp = cp[(cp >= l_boundary) & (cp <= r_boundary)]
  
  
  pi_zero = l_boundary
  pi_r_plus_one = r_boundary
  pi= c(pi_zero,I_cp,pi_r_plus_one)
  
  w = sapply(2:length(pi), function(i) integrate(f_weight,
                                                 lower = pi[i-1],
                                                 upper = pi[i])$value)
  # obtain the approximation value c
  r = length(w)
  w1 = cumsum(w)[-r]
  w2 = cumsum(w[r:1])[-r]
  ind = which(w1-w2[(r-1):1]>=0)[1]
  pi_c = ifelse(is.na(ind),
                pi_r_plus_one,
                I_cp[ind])
  c = fapprox(cp,cq,pi_c)
  
  
  ci = fapprox(cp,cq,pi)
  
  e = sum(abs(ci[1:r] - c)^l *w)
  
  return(c(e,pi_c,c))
}



####################################################

#---- 1. Evaluate the supply curve approximations----                    

####################################################


day_ahead_supply = day_ahead_supply %>%
  group_by(date,hour) %>%
  mutate(time = cur_group_id()) %>%
  ungroup() %>%
  select(p = price,
         q = quantity,
         q_cumsum = cumsum_quantity,
         time)

day_ahead_supply$p_scaled = day_ahead_supply$p/max(day_ahead_supply$p)
day_ahead_supply$q_cumsum_scaled = day_ahead_supply$q_cumsum / max(day_ahead_supply$q_cumsum)
day_ahead_supply$q_scaled =  day_ahead_supply$q / max(day_ahead_supply$q_cumsum)


##### 1.1 Approximation by L2 + marginal distribution####
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_original_supply_marginal_Approx_50nodes.RData")
n_node = 50
probabilities= seq(0, 1, length.out = n_node)
marginal_node = quantile(day_ahead_supply$p_scaled,probabilities,type = 1)
marginal_node

###### 1.1.1 L1 norm error ####
cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
#
L1_supply_EmprDis_err=foreach(i=25:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  
  p2 = marginal_node
  q2 = Approx[i,]
  
  
  p1 = p1[-1]
  p2 = p2[-1]
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 1,
                           w_type = "user_wfun",
                           user_fun = supply.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 

mean_L1_supply_EmprDis_err=mean(L1_supply_EmprDis_err)
mean_L1_supply_EmprDis_err 

###### 1.1.2 L2 norm error#####
cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

L2_supply_EmprDis_err=foreach(i=25:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  
  p2 = marginal_node
  q2 = Approx[i,]
  
  
  p1 = p1[-1]
  p2 = p2[-1]
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 2,
                           w_type = "user_wfun",
                           user_fun = supply.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 

mean_L2_supply_EmprDis_err=mean(L2_supply_EmprDis_err)
mean_L2_supply_EmprDis_err




#### 1.2 Approximation by L2 + conditional marginal distribution####
load('C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_original_supply_conditional_Approx_50nodes.RData')
Q =quantile(day_ahead_supply[(day_ahead_supply$p_scaled != 0),]$q_scaled,0.125)
# Q is the 12.5% quantile of the quantitis that are included in non-zero offers

Q = unname(Q);Q


n_node = 50
probabilities= seq(0, 1, length.out = n_node)
conditional_node = quantile(day_ahead_supply[day_ahead_supply$q_scaled > Q,]$p_scaled,probabilities,type = 1)
conditional_node
Q



###### 1.2.1 L1 norm error##### 

cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

L1_supply_CondDis_err=foreach(i=25:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  
  p2 = conditional_node
  q2 = Approx[i,]
  
  
  p1 = p1[-1]
  p2 = p2[-1]
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 1,
                           w_type = "user_wfun",
                           user_fun = supply.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 

mean_L1_supply_CondDis_err=mean(L1_supply_CondDis_err)
mean_L1_supply_CondDis_err 


###### 1.2.2 L2 norm error#####
cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

L2_supply_CondDis_err=foreach(i=25:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  
  p2 = conditional_node
  q2 = Approx[i,]
  
  
  p1 = p1[-1]
  p2 = p2[-1]
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 2,
                           w_type = "user_wfun",
                           user_fun = supply.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 

mean_L2_supply_CondDis_err=mean(L2_supply_CondDis_err)
mean_L2_supply_CondDis_err




#### 1.3 Approximation by L1 + dyadic distribution ####

##### 1.3.1 L1 norm error#####
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L1+dyadic_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L1_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L1_lst_Interval_E_50nodes.RData")
n = 50

Approx_Dyadic = L1_Approx_Dyadic_50nodes

boundaries = unname(lst_Interval_E[[n]][,2:3])
boundaries = boundaries[order(boundaries[,1]),]
boundaries[which.max(boundaries)] = Inf

unique_prices = unique(day_ahead_supply$p_scaled)
unique_boundaries = unique(c(boundaries))

real_boundaries = c()
for(b in unique_boundaries){
  pos = which.min(abs(b-unique_prices))
  real_b = unique_prices[pos]
  real_boundaries = c(real_boundaries, real_b)
}

real_boundaries[length(real_boundaries)] = Inf
boundaries[,1] = real_boundaries[1:n]
boundaries[,2] = real_boundaries[2:(n+1)]


cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
#
Approx_Dyadic_err=foreach(i=25:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
  p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  p1 = p1[-1]
  q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  
  p2 = sort(unique(c(boundaries)))
  p2 = p2[2:(length(p2)-1)]
  q2 = Approx_Dyadic[[i]][,3]
  #q2[1] = 0 # in th Spanish market case,the q of initial point is far from 0
  
  
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 1,
                           w_type = "user_wfun",
                           user_fun = supply.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 
toc()
mean_Approx_Dyadic_err_nodes=mean(Approx_Dyadic_err)
mean_Approx_Dyadic_err_nodes


##### 1.3.2 L2 norm error#####
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L1+dyadic_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L1_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L1_lst_Interval_E_50nodes.RData")
n = 50

Approx_Dyadic = L1_Approx_Dyadic_50nodes

boundaries = unname(lst_Interval_E[[n]][,2:3])
boundaries = boundaries[order(boundaries[,1]),]
boundaries[which.max(boundaries)] = Inf

unique_prices = unique(day_ahead_supply$p_scaled)
unique_boundaries = unique(c(boundaries))

real_boundaries = c()
for(b in unique_boundaries){
  pos = which.min(abs(b-unique_prices))
  real_b = unique_prices[pos]
  real_boundaries = c(real_boundaries, real_b)
}

real_boundaries[length(real_boundaries)] = Inf
boundaries[,1] = real_boundaries[1:n]
boundaries[,2] = real_boundaries[2:(n+1)]


cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

Approx_Dyadic_err=foreach(i=25:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
  p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  p1 = p1[-1]
  q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  
  p2 = sort(unique(c(boundaries)))
  p2 = p2[2:(length(p2)-1)]
  q2 = Approx_Dyadic[[i]][,3]
  #q2[1] = 0 # in th Spanish market case,the q of initial point is far from 0
  
  
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 2,
                           w_type = "user_wfun",
                           user_fun = supply.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 
toc()
mean_Approx_Dyadic_err_nodes=mean(Approx_Dyadic_err)
mean_Approx_Dyadic_err_nodes






#### 1.4 Approximation by L2 + dyadic distribution ####

##### 1.4.1 L1 norm error#####

load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L2+dyadic_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L2_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L2_lst_Interval_E_50nodes.RData")
n = 50

Approx_Dyadic = L2_Approx_Dyadic_50nodes

boundaries = unname(lst_Interval_E[[n]][,2:3])
boundaries = boundaries[order(boundaries[,1]),]
boundaries[which.max(boundaries)] = Inf

unique_prices = unique(day_ahead_supply$p_scaled)
unique_boundaries = unique(c(boundaries))

real_boundaries = c()
for(b in unique_boundaries){
  pos = which.min(abs(b-unique_prices))
  real_b = unique_prices[pos]
  real_boundaries = c(real_boundaries, real_b)
}

real_boundaries[length(real_boundaries)] = Inf
boundaries[,1] = real_boundaries[1:n]
boundaries[,2] = real_boundaries[2:(n+1)]


cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

Approx_Dyadic_err=foreach(i=25:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
  p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  p1 = p1[-1]
  q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  
  p2 = sort(unique(c(boundaries)))
  p2 = p2[2:(length(p2)-1)]
  q2 = Approx_Dyadic[[i]][,3]
  #q2[1] = 0 # in th Spanish market case,the q of initial point is far from 0
  
  
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 1,
                           w_type = "user_wfun",
                           user_fun = supply.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 
toc()
mean_Approx_Dyadic_err_nodes=mean(Approx_Dyadic_err)
mean_Approx_Dyadic_err_nodes


##### 1.4.2 L2 norm error#####
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L2+dyadic_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L2_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_supply_L2_lst_Interval_E_50nodes.RData")
n = 50

Approx_Dyadic = L2_Approx_Dyadic_50nodes

boundaries = unname(lst_Interval_E[[n]][,2:3])
boundaries = boundaries[order(boundaries[,1]),]
boundaries[which.max(boundaries)] = Inf

unique_prices = unique(day_ahead_supply$p_scaled)
unique_boundaries = unique(c(boundaries))

real_boundaries = c()
for(b in unique_boundaries){
  pos = which.min(abs(b-unique_prices))
  real_b = unique_prices[pos]
  real_boundaries = c(real_boundaries, real_b)
}

real_boundaries[length(real_boundaries)] = Inf
boundaries[,1] = real_boundaries[1:n]
boundaries[,2] = real_boundaries[2:(n+1)]


cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

Approx_Dyadic_err=foreach(i=25:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
  p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  p1 = p1[-1]
  q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  
  p2 = sort(unique(c(boundaries)))
  p2 = p2[2:(length(p2)-1)]
  q2 = Approx_Dyadic[[i]][,3]
  #q2[1] = 0 # in th Spanish market case,the q of initial point is far from 0
  
  
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 2,
                           w_type = "user_wfun",
                           user_fun = supply.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 
toc()
mean_Approx_Dyadic_err_nodes=mean(Approx_Dyadic_err)
mean_Approx_Dyadic_err_nodes
















####################################################

#---- 2. Evaluate the demand curve approximations----                    

####################################################

day_ahead_demand = day_ahead_demand %>%
  group_by(date,hour) %>%
  mutate(time = cur_group_id()) %>%
  ungroup() %>%
  select(p = price,
         q = quantity,
         q_cumsum = cumsum_quantity,
         time)

day_ahead_demand$p_scaled = day_ahead_demand$p/max(day_ahead_demand$p)
day_ahead_demand$q_cumsum_scaled = day_ahead_demand$q_cumsum / max(day_ahead_demand$q_cumsum)
day_ahead_demand$q_scaled =  day_ahead_demand$q / max(day_ahead_demand$q_cumsum) 
# must be max(day_ahead_demand$q_cumsum), otherwise, curves end at differnt height




##### 2.1 Approximation by L2 + marginal distribution####

load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_original_demand_marginal_Approx_50nodes.RData")

n_node = 50
probabilities= seq(0, 1, length.out = n_node)
marginal_node = quantile(day_ahead_demand$p_scaled,probabilities,type = 1)
marginal_node


###### 2.1.1 L1 norm error ####


cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

L1_demand_EmprDis_err=foreach(i=25:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  p2 = marginal_node
  q2 = Approx[i,]
  #q2[1] = 0 # in th Spanish market case,the q of initial point is far from 0
  
  p1 = p1[-1]
  p2 = p2[-1]
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 1,
                           w_type = "user_wfun",
                           user_fun = demand.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 

mean_L1_demand_EmprDis_err=mean(L1_demand_EmprDis_err)
mean_L1_demand_EmprDis_err 




###### 2.1.2 L2 norm error#####


cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

L2_demand_EmprDis_err=foreach(i=25:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  p2 = marginal_node
  q2 = Approx[i,]
  #q2[1] = 0 # in th Spanish market case,the q of initial point is far from 0
  
  p1 = p1[-1]
  p2 = p2[-1]
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 2,
                           w_type = "user_wfun",
                           user_fun = demand.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 

mean_L2_demand_EmprDis_err=mean(L2_demand_EmprDis_err)
mean_L2_demand_EmprDis_err 





#### 2.2 Approximation by L2 + conditional marginal distribution####
load('C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_original_demand_conditional_Approx_50nodes.RData')
Q =quantile(day_ahead_demand[(day_ahead_demand$p_scaled != 0),]$q_scaled,0.125)
# Q is the 12.5% quantile of the quantitis that are included in non-zero offers

Q = unname(Q);Q

n_node = 50
probabilities= seq(0, 1, length.out = n_node)
conditional_node = quantile(day_ahead_demand[day_ahead_demand$q_scaled > Q,]$p_scaled,probabilities,type = 1)
conditional_node
Q



###### 2.2.1 L1 norm error##### 

cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

L1_demand_CondDis_err=foreach(i=25:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  p2 = conditional_node
  q2 = Approx[i,]
  
  
  p1 = p1[-1]
  p2 = p2[-1]
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 1,
                           w_type = "user_wfun",
                           user_fun = demand.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 

mean_L1_demand_CondDis_err=mean(L1_demand_CondDis_err)
mean_L1_demand_CondDis_err 


###### 2.2.2 L2 norm error#####
cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

L2_demand_CondDis_err=foreach(i=25:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  p2 = conditional_node
  q2 = Approx[i,]
  
  
  p1 = p1[-1]
  p2 = p2[-1]
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 2,
                           w_type = "user_wfun",
                           user_fun = demand.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 

mean_L2_demand_CondDis_err=mean(L2_demand_CondDis_err)
mean_L2_demand_CondDis_err


#### 2.3 Approximation by L1 + dyadic distribution ####

##### 2.3.1 L1 norm error#####
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L1+dyadic_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L1_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L1_lst_Interval_E_50nodes.RData")
n = 50

Approx_Dyadic = L1_Approx_Dyadic_50nodes

boundaries = unname(lst_Interval_E[[n]][,2:3])
boundaries = boundaries[order(boundaries[,1]),]
boundaries[which.max(boundaries)] = Inf

unique_prices = unique(day_ahead_demand$p_scaled)
unique_boundaries = unique(c(boundaries))

real_boundaries = c()
for(b in unique_boundaries){
  pos = which.min(abs(b-unique_prices))
  real_b = unique_prices[pos]
  real_boundaries = c(real_boundaries, real_b)
}

real_boundaries[length(real_boundaries)] = Inf
boundaries[,1] = real_boundaries[1:n]
boundaries[,2] = real_boundaries[2:(n+1)]


cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

Approx_Dyadic_err=foreach(i=25:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
  p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  p1 = p1[-1]
  q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  p2 = sort(unique(c(boundaries)))
  p2 = p2[2:(length(p2)-1)]
  q2 = Approx_Dyadic[[i]][,3]
  #q2[1] = 0 # in th Spanish market case,the q of initial point is far from 0
  
  
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 1,
                           w_type = "user_wfun",
                           user_fun = demand.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 
toc()
mean_Approx_Dyadic_err_nodes=mean(Approx_Dyadic_err)
mean_Approx_Dyadic_err_nodes


##### 2.3.2 L2 norm error#####
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L1+dyadic_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L1_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L1_lst_Interval_E_50nodes.RData")
n = 50

Approx_Dyadic = L1_Approx_Dyadic_50nodes

boundaries = unname(lst_Interval_E[[n]][,2:3])
boundaries = boundaries[order(boundaries[,1]),]
boundaries[which.max(boundaries)] = Inf

unique_prices = unique(day_ahead_demand$p_scaled)
unique_boundaries = unique(c(boundaries))

real_boundaries = c()
for(b in unique_boundaries){
  pos = which.min(abs(b-unique_prices))
  real_b = unique_prices[pos]
  real_boundaries = c(real_boundaries, real_b)
}

real_boundaries[length(real_boundaries)] = Inf
boundaries[,1] = real_boundaries[1:n]
boundaries[,2] = real_boundaries[2:(n+1)]


cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

Approx_Dyadic_err=foreach(i=25:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
  p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  p1 = p1[-1]
  q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  p2 = sort(unique(c(boundaries)))
  p2 = p2[2:(length(p2)-1)]
  q2 = Approx_Dyadic[[i]][,3]
  #q2[1] = 0 # in th Spanish market case,the q of initial point is far from 0
  
  
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 2,
                           w_type = "user_wfun",
                           user_fun = demand.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 
toc()
mean_Approx_Dyadic_err_nodes=mean(Approx_Dyadic_err)
mean_Approx_Dyadic_err_nodes






#### 2.4 Approximation by L2 + dyadic distribution ####

##### 2.4.1 L1 norm error#####

load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L2+dyadic_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L2_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L2_lst_Interval_E_50nodes.RData")
n = 50

Approx_Dyadic = L2_Approx_Dyadic_50nodes

boundaries = unname(lst_Interval_E[[n]][,2:3])
boundaries = boundaries[order(boundaries[,1]),]
boundaries[which.max(boundaries)] = Inf

unique_prices = unique(day_ahead_demand$p_scaled)
unique_boundaries = unique(c(boundaries))

real_boundaries = c()
for(b in unique_boundaries){
  pos = which.min(abs(b-unique_prices))
  real_b = unique_prices[pos]
  real_boundaries = c(real_boundaries, real_b)
}

real_boundaries[length(real_boundaries)] = Inf
boundaries[,1] = real_boundaries[1:n]
boundaries[,2] = real_boundaries[2:(n+1)]


cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

Approx_Dyadic_err=foreach(i=25:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
  p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  p1 = p1[-1]
  q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  p2 = boundaries[,1]
  p2 = p2[-1]
  q2 = Approx_Dyadic[[i]][,3]
  #q2[1] = 0 # in th Spanish market case,the q of initial point is far from 0
  
  
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 1,
                           w_type = "user_wfun",
                           user_fun = demand.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 
toc()
mean_Approx_Dyadic_err_nodes=mean(Approx_Dyadic_err)
mean_Approx_Dyadic_err_nodes


##### 2.4.2 L2 norm error#####
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L2+dyadic_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L2_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/sp_day-ahead_demand_L2_lst_Interval_E_50nodes.RData")
n = 50

Approx_Dyadic = L2_Approx_Dyadic_50nodes

boundaries = unname(lst_Interval_E[[n]][,2:3])
boundaries = boundaries[order(boundaries[,1]),]
boundaries[which.max(boundaries)] = Inf

unique_prices = unique(day_ahead_demand$p_scaled)
unique_boundaries = unique(c(boundaries))

real_boundaries = c()
for(b in unique_boundaries){
  pos = which.min(abs(b-unique_prices))
  real_b = unique_prices[pos]
  real_boundaries = c(real_boundaries, real_b)
}

real_boundaries[length(real_boundaries)] = Inf
boundaries[,1] = real_boundaries[1:n]
boundaries[,2] = real_boundaries[2:(n+1)]


cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

Approx_Dyadic_err=foreach(i=25:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
  p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  p1 = p1[-1]
  q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  p2 = boundaries[,1]
  p2 = p2[-1]
  q2 = Approx_Dyadic[[i]][,3]
  #q2[1] = 0 # in th Spanish market case,the q of initial point is far from 0
  
  
  stepwise_distance_w_plot(p1 = p1, 
                           q1 = q1,
                           p2 = p2,
                           q2 = q2,
                           l = 2,
                           w_type = "user_wfun",
                           user_fun = demand.mixture.gaussian,
                           plot = F)
}
close(pb)
stopCluster(cl) 
toc()
mean_Approx_Dyadic_err_nodes=mean(Approx_Dyadic_err)
mean_Approx_Dyadic_err_nodes



