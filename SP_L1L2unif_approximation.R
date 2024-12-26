library(tidyverse)
library(cubature)
library(akima)
library(nor1mix)
library(LaplacesDemon)
library(mixtools)

library(foreach)
library(doParallel)
library(doSNOW)

load("C:/Zehang_workspace/Approximation_Supply_Demand/day_ahead_original_supply_curves.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/day_ahead_original_demand_curves.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gammamix_match_prices.RData")

#############################################################################
#                     0. Defining some useful functions                     #
#############################################################################

#### 0.1 distance function ####
stepwise_distance_w_plot <- function(p1, q1, #prices and quantities for curve 1
                                     p2, q2, #prices and quantities for curve 2
                                     l , #l of the distance
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

#### 0.2 weight function being used in the SUPPLY CURVE approximation procedure ####
dgammamixEM <- function(x,GammaDist) {
  dist = sapply(1:length(GammaDist$lambda),function(i) GammaDist$lambda[i]*dgamma(x, GammaDist$gamma.pars[1,i], 
                                                                                  1/GammaDist$gamma.pars[2,i]))
  return(sum(dist))
}

supply.mixture.gamma =  function(x) dgammamixEM(x,gammamixEM_p_m2)

#supply.mixture.gaussian = function(x) ifelse(abs(x-1/2)<=1/2,1,0)

# load("C:/Zehang_workspace/Approximation_Supply_Demand/norMixmodel_supply_p_scaled.RData")
# supply.mixture.gaussian = function(x) dnorMix(norMixmodel_supply_p_scaled,x)


# supply_p_scaled = day_ahead_supply$price/max(day_ahead_supply$price)
# norMixmodel_supply_p_scaled <-  norMixMLE(supply_p_scaled,m=3)
# plot(density(supply_p_scaled),col=3)
# lines(norMixmodel_supply_p_scaled)


#### 0.3 weight function being used in the DEMAND CURVE approximation procedure ####
dgammamixEM <- function(x,GammaDist) {
  dist = sapply(1:length(GammaDist$lambda),function(i) GammaDist$lambda[i]*dgamma(x, GammaDist$gamma.pars[1,i], 
                                                                                  1/GammaDist$gamma.pars[2,i]))
  return(sum(dist))
}

demand.mixture.gamma =  function(x) dgammamixEM(x,gammamixEM_p_m2)


#demand.mixture.gaussian = function(x) ifelse(abs(x-1/2)<=1/2,1,0)


# load("C:/Zehang_workspace/Approximation_Supply_Demand/norMixmodel_demand_p_scaled.RData")
# demand.mixture.gaussian = function(x) dnorMix(norMixmodel_demand_p_scaled,x)


# demand_p_scaled = day_ahead_demand$price / max(day_ahead_demand$price)
# norMixmodel_demand_p_scaled = norMixMLE(demand_p_scaled,m=3)
# plot(density(demand_p_scaled),col=3)
# lines(norMixmodel_demand_p_scaled)


#### 0.4 functions calculating CW, W and approximation ####
f_CW = function(p,q,lower,upper,f_weight) {
  p = p[-1]
  #f_Cp = stepfun(p, q, f = 1, right = T)
  f_Cp = stepfun(p, q, f = 0)
  CW = function(x) f_Cp(x) * f_weight(x)
  return(cubintegrate(CW, 
                      upper = upper, 
                      lower = lower,method = "pcubature")$integral)
}

f_W = function(lower,upper,f_weight) cubintegrate(f_weight, upper = upper, lower = lower,method = "pcubature")$integral

# Approximating function
f_L2Approx = function(i,c_vector,node,p){c_vector[i]* (p >= node[i])}

#### 0.5 functions calculating approximation and setting nodes####
fapprox = function(x,y,eval) {
  x = x[-1]
  sf = stepfun(x, y, f = 0) 
  return(sf(eval))}

f_interval_approx_error = function(cp,cq, # the p and q of all jumps
                                   l_boundary, # the left boundary of interval
                                   r_boundary, # the right boundary of interval
                                   f_weight,
                                   increasing, # increasing = T for supply curves, F for demand curves
                                   l ){# L1 or L2 norm
  
  I_cp = cp[(cp > l_boundary) & (cp < r_boundary)] #I_cp = cp[(cp >= l_boundary) & (cp < r_boundary)]
  
  
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
  ind = ifelse(increasing == T,which(w1-w2[(r-1):1]>=0)[1],tail(which(-w1+w2[(r-1):1]>=0),1))
  pi_c = ifelse(is.na(ind),
                ifelse(increasing == T,pi_r_plus_one,pi_zero), #pi_r_plus_one
                I_cp[ind]) 
  c = fapprox(cp,cq,pi_c)
  
  
  ci = fapprox(cp,cq,pi)
  
  e = sum(abs(ci[1:r] - c)^l *w)
  
  return(c(e,pi_c,c))
}


E = matrix(NA,nrow = 0,ncol=5)
colnames(E) = c('error','approx_p','approx_q','l_boundary','r_boundary')


max_q_cumsum = 90834.4#max(max(day_ahead_demand$q_cumsum),max(day_ahead_supply$q_cumsum))
max_p = 180.3#max(max(day_ahead_demand$p),max(day_ahead_supply$p))

#############################################################################
#                 1. L2 Approximation of supply curves                      #
#############################################################################

day_ahead_supply = day_ahead_supply %>%
  group_by(date,hour) %>%
  mutate(time = cur_group_id()) %>%
  ungroup() %>%
  dplyr::select(p = price,
                q = quantity,
                q_cumsum = cumsum_quantity,
                time)

day_ahead_supply$p_scaled = day_ahead_supply$p/max_p
day_ahead_supply$q_cumsum_scaled = day_ahead_supply$q_cumsum / max_q_cumsum
day_ahead_supply$q_scaled =  day_ahead_supply$q / max_q_cumsum
# must be max(day_ahead_supply$q_cumsum), otherwise, curves end at differnt height

n_node = 50
node= seq(0, 1, length.out = n_node)
node

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
#
Approx = foreach(i=1:max(day_ahead_supply$time),.combine='rbind', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  
  numerator = sapply(1:(length(node)-1), function(n) f_CW(p1,q1,lower = node[-length(node)][n],upper = node[-1][n],f_weight=supply.mixture.gamma))
  denominator = sapply(1:(length(node)-1), function(n) f_W(lower = node[-length(node)][n],upper = node[-1][n],f_weight = supply.mixture.gamma))
  
  frac = numerator/denominator
  cn = f_CW(p1,q1,lower = node[length(node)-1], upper=node[length(node)], supply.mixture.gamma)/ f_W(lower = node[length(node)-1], upper=node[length(node)],supply.mixture.gamma)
  
  c_vector = c(frac, cn) - c(0,frac)           
  
  Approx = sapply(1:length(node), function(n) sapply(node, function(p) f_L2Approx(n,c_vector,node,p)))
  Approx = rowSums(Approx)
  Approx
}
close(pb)
stopCluster(cl) 

save_path = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_L2_Unif_GammaWeight_Approx_",n_node,"nodes",".RData",sep = "")
save(Approx,file=save_path)

i = sample(1:43848,1)
p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
p1 = p1[-1]
q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled

p2 = node
q2 = Approx[i,]
p2 = p2[-1]

sfun1 = stepfun(p1, q1, f = 0)
sfun2 = stepfun(p2, q2, f = 0)

plot(sfun1,xlim = c(0,1.1),
     ylab = "Quantity (MW)", xlab = paste('Price (', "\u20AC", '/MW)',sep=''),
     main = '',xaxs="i",yaxs="i", pch = 16)
lines(sfun2,col='blue',xlim = c(0,1.1),pch = 16)

#############################################################################
#                 2. L1 Approximation of supply curves                      #
#############################################################################


n_node = 50
node= seq(0, 1, length.out = n_node+1)

boundaries = cbind(node[1:50],node[2:51])

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
#1:max(day_ahead_supply$time)
L1_Approx_Unif_50nodes = foreach(i=1:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
  cp = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  cq = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  
  res = sapply(1:n_node, function(m) f_interval_approx_error(cp,cq,
                                                             boundaries[m,1],boundaries[m,2], 
                                                             Vectorize(supply.mixture.gamma),increasing = T,
                                                             l=1))
  list(t(res))
  
}
close(pb)
stopCluster(cl) 

save_path = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_L1_Unif_GammaWeight_Approx_",n_node,"nodes",".RData",sep = "")
save(L1_Approx_Unif_50nodes,file=save_path)

i =  sample(1:43848,1)
p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
p1 = p1[-1]
q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled

p2 = sort(unique(c(boundaries)))
p2 = p2[2:(length(p2)-1)]
q2 = L1_Approx_Unif_50nodes[[i]][,3]

sfun1 = stepfun(p1, q1, f = 0)
sfun2 = stepfun(p2, q2, f = 0)
plot(sfun1,xlim = c(0,1.1), ylim = c(0,1),
     ylab = "Quantity (MW)", xlab = paste('Price (', "\u20AC", '/MW)',sep=''),
     main = '',xaxs="i",yaxs="i", pch = 16)
lines(sfun2,col='blue',xlim = c(0,1.1),pch = 16)



#############################################################################
#                    3. Approximation of demand curves                      #
#############################################################################


day_ahead_demand = day_ahead_demand %>%
  group_by(date,hour) %>%
  mutate(time = cur_group_id()) %>%
  ungroup() %>%
  dplyr::select(p = price,
         q = quantity,
         q_cumsum = cumsum_quantity,
         time)

day_ahead_demand$p_scaled = day_ahead_demand$p/max_p
day_ahead_demand$q_cumsum_scaled = day_ahead_demand$q_cumsum / max_q_cumsum
day_ahead_demand$q_scaled =  day_ahead_demand$q / max_q_cumsum


n_node = 50
node= seq(0, 1, length.out = n_node)
node

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
#
Approx = foreach(i=1:max(day_ahead_demand$time),.combine='rbind', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  Approx = fapprox(p1,q1,node)
  Approx
}
close(pb)
stopCluster(cl) 

save_path = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_L2_Unif_GammaWeight_Approx_",n_node,"nodes",".RData",sep = "")
save(Approx,file=save_path)

i = sample(1:43848,1)
p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
p1 = p1[-1]
q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled

p2 = node
q2 = Approx[i,]
p2 = p2[-1]

sfun1 = stepfun(p1, q1, f = 0)
sfun2 = stepfun(p2, q2, f = 0)

plot(sfun1,xlim = c(0,1.1),
     ylab = "Quantity (MW)", xlab = paste('Price (', "\u20AC", '/MW)',sep=''),
     main = '',xaxs="i",yaxs="i", pch = 16)
lines(sfun2,col='blue',xlim = c(0,1.1),pch = 16)



#############################################################################
#                 4. L1 Approximation of demand curves                      #
#############################################################################


n_node = 50
node= seq(0, 1, length.out = n_node+1)

boundaries = cbind(node[1:50],node[2:51])

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
#
L1_Approx_Unif_50nodes = foreach(i=1:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
  cp = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  cq = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  res = sapply(1:n_node, function(m) f_interval_approx_error(cp,cq,
                                                             boundaries[m,1],boundaries[m,2], 
                                                             Vectorize(demand.mixture.gamma),increasing = F,
                                                             l=1))
  list(t(res))
  
}
close(pb)
stopCluster(cl) 

save_path = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_L1_Unif_GammaWeight_Approx_",n_node,"nodes",".RData",sep = "")
save(L1_Approx_Unif_50nodes,file=save_path)

i =  sample(1:43848,1)
p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
p1 = p1[-1]
q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled

p2 = sort(unique(c(boundaries)))
p2 = p2[2:(length(p2)-1)]
q2 = L1_Approx_Unif_50nodes[[i]][,3]

sfun1 = stepfun(p1, q1, f = 0)
sfun2 = stepfun(p2, q2, f = 0)
plot(sfun1,xlim = c(0,1.1), 
     ylab = "Quantity (MW)", xlab = paste('Price (', "\u20AC", '/MW)',sep=''),
     main = '',xaxs="i",yaxs="i", pch = 16)
lines(sfun2,col='blue',xlim = c(0,1.1),pch = 16)

