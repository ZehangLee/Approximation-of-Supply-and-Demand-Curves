library(tidyverse)
library(cubature)
library(akima)
library(nor1mix)
library(LaplacesDemon)

library(foreach)
library(doParallel)
library(doSNOW)
library(tictoc)

load("C:/Zehang_workspace/Approximation_Supply_Demand/day_ahead_original_supply_curves.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/day_ahead_original_demand_curves.RData")
load('C:/Zehang_workspace/Approximation_Supply_Demand/gammamix_match_prices.RData')


#############################################################################
#                     0. Defining some useful functions                     #
#############################################################################


#### 0.1 weight function being used in the SUPPLY CURVE approximation procedure ####
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
# norMixmodel_supply_p_scaled <-  norMixMLE(supply_p_scaled,m=4)
# plot(density(supply_p_scaled),col=3)
# lines(norMixmodel_supply_p_scaled)


#### 0.2 weight function being used in the DEMAND CURVE approximation procedure ####
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
# norMixmodel_demand_p_scaled = norMixMLE(demand_p_scaled,m=4)
# plot(density(demand_p_scaled),col=3)
# lines(norMixmodel_demand_p_scaled)


#### 0.3 functions calculating approximation and setting nodes####
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
  
  I_cp = cp[(cp > l_boundary) & (cp < r_boundary)]
  
  
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
                ifelse(increasing == T,pi_r_plus_one,pi_zero),
                I_cp[ind]) 
  c = fapprox(cp,cq,pi_c)
  
  
  ci = fapprox(cp,cq,pi)
  
  e = sum(abs(ci[1:r] - c)^l *w)
  
  return(c(e,pi_c,c))
}


E = matrix(NA,nrow = 0,ncol=5)
colnames(E) = c('error','approx_p','approx_q','l_boundary','r_boundary')


#### 0.4 functions calculating L2 approximation####
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

max_q_cumsum = 90834.4#max(max(day_ahead_demand$q_cumsum),max(day_ahead_supply$q_cumsum))
max_p = 180.3


#############################################################################
#                    1. Approximation of supply curves                      #
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

#### 1.1 Marginal + L1 dyadic ####

#### 1.1.1 Prior (n_nodes)/2 nodes placed by distribution and the remainder by dyadic ####

n_nodes = 50
n_prior_node = n_nodes/2 
probabilities= seq(0, 1, length.out = n_prior_node+1)[1:(n_prior_node)]
prior_node = quantile(day_ahead_supply$p_scaled,probabilities,type = 1)
prior_node = unname(prior_node)
prior_node

prior_l_boundaries = prior_node
prior_r_boundaries = c(prior_node[2:length(prior_node)],1)
prior_boundaries = cbind(prior_l_boundaries,prior_r_boundaries)

n_later_node = n_nodes/2


n_curve = max(day_ahead_supply$time)

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = n_curve
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

tic()
prior_error = foreach(cc=1:n_curve,.options.snow=opts)%dopar%{
  
  cp_real = day_ahead_supply[day_ahead_supply$time==cc,]$p_scaled
  cq_real = day_ahead_supply[day_ahead_supply$time==cc,]$q_cumsum_scaled
  Approx = sapply(1:n_prior_node, function(m) f_interval_approx_error(cp_real,cq_real,
                                                                      prior_boundaries[m,1],prior_boundaries[m,2], 
                                                                      Vectorize(supply.mixture.gamma),
                                                                      increasing = T,
                                                                      l=1))
  
}

close(pb)
stopCluster(cl)

Interval_E = cbind(Reduce("+", prior_error)[1,]/n_curve,
                   prior_l_boundaries,
                   prior_r_boundaries)
colnames(Interval_E) = c('error','l_boundary','r_boundary')
lst_CurvE = lapply(prior_error, function(x) {x = t(x)
x = cbind(x,prior_l_boundaries,prior_r_boundaries)
colnames(x) = c('error','approx_p','approx_q','l_boundary','r_boundary')
x})
lst_Interval_E = list()


for(n in 1:n_later_node){
  
  p1 = unname(Interval_E[which.max(Interval_E[,1]),2])
  p3 = unname(Interval_E[which.max(Interval_E[,1]),3])
  p2 = unname((p1+p3)/2)
  
  Interval_drop = which.max(Interval_E[,1])
  Interval_E = Interval_E[-Interval_drop,]
  
  for(cc in 1:n_curve){
    
    print(paste0('n_later_node=',n,'; ','curve=',cc))
    E= lst_CurvE[[cc]]
    E = E[-Interval_drop,]
    
    cp = day_ahead_supply[day_ahead_supply$time==cc,]$p_scaled
    cq = day_ahead_supply[day_ahead_supply$time==cc,]$q_cumsum_scaled
    
    Il = c(p1,p2)
    Ir = c(p2,ifelse(p3 == max(cp),Inf,p3))
    
    
    left = f_interval_approx_error(cp,cq,Il[1],Il[2], Vectorize(supply.mixture.gamma),l=1,increasing = T)
    right = f_interval_approx_error(cp,cq,Ir[1],Ir[2], Vectorize(supply.mixture.gamma),l=1,increasing = T)
    
    
    E = rbind(E,
              c(left,p1,p2),
              c(right,p2,p3))
    
    lst_CurvE[[cc]] = E
    
  }
  
  
  Interval_mean_error = Reduce("+", lst_CurvE)[,1]/n_curve
  Interval_mean_error = tail(Interval_mean_error,2)
  Interval_E = rbind(Interval_E,
                     c(Interval_mean_error[1],p1,p2),
                     c(Interval_mean_error[2],p2,p3))
  lst_Interval_E[[n]] = Interval_E
  
}
toc()

save_path1 = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L1_lst_CurvE_",n_nodes,"nodes",".RData",sep = "")
save(lst_CurvE,file=save_path1)
save_path2 = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L1_lst_Interval_E_",n_nodes,"nodes",".RData",sep = "")
save(lst_Interval_E,file = save_path2)

##### 1.1.2. obtain the approximation with 50 nodes ####
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L1_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L1_lst_Interval_E_50nodes.RData")
n_node = 50
n = n_node/2

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
boundaries[,1] = real_boundaries[1:n_node]
boundaries[,2] = real_boundaries[2:(n_node+1)]


cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

tic()
L1_Approx_Dyadic_50nodes=foreach(i=1:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
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
toc()

save_path = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L1+dyadic_",n,"nodes",".RData",sep = "")
save(L1_Approx_Dyadic_50nodes,file=save_path)


i =  sample(1:43848,1)
p1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
p1 = p1[-1]
q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled

p2 = sort(unique(c(boundaries)))
p2 = p2[2:(length(p2)-1)]
q2 = L1_Approx_Dyadic_50nodes[[i]][,3]

sfun1 = stepfun(p1, q1, f = 0)
sfun2 = stepfun(p2, q2, f = 0)
plot(sfun1,xlim = c(0,1.1), ylim = c(0,1),
     ylab = "Quantity (MW)", xlab = paste('Price (', "\u20AC", '/MW)',sep=''),
     main = '',xaxs="i",yaxs="i", pch = 16)
lines(sfun2,col='blue',xlim = c(0,1.1),pch = 16)


#### 1.2 Marginal + L2 + dyadic ####

#### 1.2.1 Prior (n_nodes)/2 nodes placed by distribution and the remainder by dyadic ####
n_nodes = 50
n_prior_node = n_nodes/2 
probabilities= seq(0, 1, length.out = n_prior_node+1)[1:(n_prior_node)]
prior_node = quantile(day_ahead_supply$p_scaled,probabilities,type = 1)
prior_node = unname(prior_node)
prior_node

prior_l_boundaries = prior_node
prior_r_boundaries = c(prior_node[2:length(prior_node)],1)
prior_boundaries = cbind(prior_l_boundaries,prior_r_boundaries)

n_later_node = n_nodes/2


n_curve = max(day_ahead_supply$time)

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = n_curve
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

tic()
prior_error = foreach(cc=1:n_curve,.options.snow=opts)%dopar%{
  
  cp_real = day_ahead_supply[day_ahead_supply$time==cc,]$p_scaled
  cq_real = day_ahead_supply[day_ahead_supply$time==cc,]$q_cumsum_scaled
  Approx = sapply(1:n_prior_node, function(m) f_interval_approx_error(cp_real,cq_real,
                                                                      prior_boundaries[m,1],prior_boundaries[m,2], 
                                                                      Vectorize(supply.mixture.gamma),
                                                                      increasing = T,
                                                                      l=2))
  
}

close(pb)
stopCluster(cl)

Interval_E = cbind(Reduce("+", prior_error)[1,]/n_curve,
                   prior_l_boundaries,
                   prior_r_boundaries)
colnames(Interval_E) = c('error','l_boundary','r_boundary')
lst_CurvE = lapply(prior_error, function(x) {x = t(x)
x = cbind(x,prior_l_boundaries,prior_r_boundaries)
colnames(x) = c('error','approx_p','approx_q','l_boundary','r_boundary')
x})
lst_Interval_E = list()

#
for(n in 1:n_later_node){
  
  p1 = unname(Interval_E[which.max(Interval_E[,1]),2])
  p3 = unname(Interval_E[which.max(Interval_E[,1]),3])
  p2 = unname((p1+p3)/2)
  
  Interval_drop = which.max(Interval_E[,1])
  Interval_E = Interval_E[-Interval_drop,]
  #n_curve
  for(cc in 1:n_curve){
    
    print(paste0('n_later_node=',n,'; ','curve=',cc))
    E= lst_CurvE[[cc]]
    E = E[-Interval_drop,]
    
    cp = day_ahead_supply[day_ahead_supply$time==cc,]$p_scaled
    cq = day_ahead_supply[day_ahead_supply$time==cc,]$q_cumsum_scaled
    
    Il = c(p1,p2)
    Ir = c(p2,ifelse(p3 == max(cp),Inf,p3))
    
    
    left = f_interval_approx_error(cp,cq,Il[1],Il[2], Vectorize(supply.mixture.gamma),l=2,increasing = T)
    right = f_interval_approx_error(cp,cq,Ir[1],Ir[2], Vectorize(supply.mixture.gamma),l=2,increasing = T)
    
    
    E = rbind(E,
              c(left,p1,p2),
              c(right,p2,p3))
    
    lst_CurvE[[cc]] = E
    
  }
  
  
  Interval_mean_error = Reduce("+", lst_CurvE)[,1]/n_curve
  Interval_mean_error = tail(Interval_mean_error,2)
  Interval_E = rbind(Interval_E,
                     c(Interval_mean_error[1],p1,p2),
                     c(Interval_mean_error[2],p2,p3))
  lst_Interval_E[[n]] = Interval_E
  
}

toc()

save_path1 = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L2_lst_CurvE_",n_nodes,"nodes",".RData",sep = "")
save(lst_CurvE,file=save_path1)
save_path2 = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L2_lst_Interval_E_",n_nodes,"nodes",".RData",sep = "")
save(lst_Interval_E,file = save_path2)

##### 1.2.2. obtain the approximation with 50 nodes ####
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L2_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L2_lst_Interval_E_50nodes.RData")

n_node = 50
n = n_node/2

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
boundaries[,1] = real_boundaries[1:n_node]
boundaries[,2] = real_boundaries[2:(n_node+1)]

node = boundaries[,1]


cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

tic()
L2_Approx_Dyadic_50nodes = foreach(i=1:max(day_ahead_supply$time),.combine='rbind', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
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
toc()

save_path = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L2+dyadic_",n,"nodes2",".RData",sep = "")
save(L2_Approx_Dyadic_50nodes,file=save_path)

i = sample(1:43848,1)
cp1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
cq1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled

cp2 = node
cq2 = L2_Approx_Dyadic_50nodes[i,]

cp1 = cp1[-1]
cp2 = cp2[-1]
sfun1 = stepfun(cp1, cq1, f = 0)
sfun2 = stepfun(cp2, cq2, f = 0)

plot(sfun1, xlim = c(0, max(cp1,cp2)), ylim = c(min(cq1, cq2), max(cq1, cq2)),
     ylab = "Scaled Quantity", xlab = "Scaled Price", 
     main = '',xaxs="i",yaxs="i",pch=16)
plot(sfun2,  col = "blue",  xlim = c(0, max(cp1,cp2)), add = TRUE,pch=16)





#############################################################################
#                    2. Approximation of demand curves                      #
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

#### 2.1 Marginal + L1 dyadic ####

#### 2.1.1 Prior (n_nodes)/2 nodes placed by distribution and the remainder by dyadic ####

n_nodes = 50
n_prior_node = n_nodes/2 
probabilities= seq(0, 1, length.out = n_prior_node+1)[1:(n_prior_node)]
prior_node = quantile(day_ahead_demand$p_scaled,probabilities,type = 1)
prior_node = unname(prior_node)
prior_node

prior_l_boundaries = prior_node
prior_r_boundaries = c(prior_node[2:length(prior_node)],1)
prior_boundaries = cbind(prior_l_boundaries,prior_r_boundaries)

n_later_node = n_nodes/2


n_curve = max(day_ahead_demand$time)

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = n_curve
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

tic()
prior_error = foreach(cc=1:n_curve,.options.snow=opts)%dopar%{
  
  cp_real = day_ahead_demand[day_ahead_demand$time==cc,]$p_scaled
  cq_real = day_ahead_demand[day_ahead_demand$time==cc,]$q_cumsum_scaled
  Approx = sapply(1:n_prior_node, function(m) f_interval_approx_error(cp_real,cq_real,
                                                                      prior_boundaries[m,1],prior_boundaries[m,2], 
                                                                      Vectorize(demand.mixture.gamma),
                                                                      increasing = F,
                                                                      l=1))
  
}

close(pb)
stopCluster(cl)

Interval_E = cbind(Reduce("+", prior_error)[1,]/n_curve,
                   prior_l_boundaries,
                   prior_r_boundaries)
colnames(Interval_E) = c('error','l_boundary','r_boundary')
lst_CurvE = lapply(prior_error, function(x) {x = t(x)
x = cbind(x,prior_l_boundaries,prior_r_boundaries)
colnames(x) = c('error','approx_p','approx_q','l_boundary','r_boundary')
x})
lst_Interval_E = list()


for(n in 1:n_later_node){
  
  p1 = unname(Interval_E[which.max(Interval_E[,1]),2])
  p3 = unname(Interval_E[which.max(Interval_E[,1]),3])
  p2 = unname((p1+p3)/2)
  
  Interval_drop = which.max(Interval_E[,1])
  Interval_E = Interval_E[-Interval_drop,]
  
  for(cc in 1:n_curve){
    
    print(paste0('n_later_node=',n,'; ','curve=',cc))
    E= lst_CurvE[[cc]]
    E = E[-Interval_drop,]
    
    cp = day_ahead_demand[day_ahead_demand$time==cc,]$p_scaled
    cq = day_ahead_demand[day_ahead_demand$time==cc,]$q_cumsum_scaled
    
    Il = c(p1,p2)
    Ir = c(p2,ifelse(p3 == max(cp),Inf,p3))
    
    
    left = f_interval_approx_error(cp,cq,Il[1],Il[2], Vectorize(demand.mixture.gamma),l=1,increasing = F)
    right = f_interval_approx_error(cp,cq,Ir[1],Ir[2], Vectorize(demand.mixture.gamma),l=1,increasing = F)
    
    
    E = rbind(E,
              c(left,p1,p2),
              c(right,p2,p3))
    
    lst_CurvE[[cc]] = E
    
  }
  
  
  Interval_mean_error = Reduce("+", lst_CurvE)[,1]/n_curve
  Interval_mean_error = tail(Interval_mean_error,2)
  Interval_E = rbind(Interval_E,
                     c(Interval_mean_error[1],p1,p2),
                     c(Interval_mean_error[2],p2,p3))
  lst_Interval_E[[n]] = Interval_E
  
}

toc()

save_path1 = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L1_lst_CurvE_",n_nodes,"nodes",".RData",sep = "")
save(lst_CurvE,file=save_path1)
save_path2 = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L1_lst_Interval_E_",n_nodes,"nodes",".RData",sep = "")
save(lst_Interval_E,file = save_path2)

##### 2.1.2. obtain the approximation with 50 nodes ####

load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L1_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L1_lst_Interval_E_50nodes.RData")

n_nodes = 50
n = n_nodes/2
n_curve = max(day_ahead_demand$time)
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
boundaries[,1] = real_boundaries[1:n_nodes]
boundaries[,2] = real_boundaries[2:(n_nodes+1)]


cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

tic()
L1_Approx_Dyadic_50nodes=foreach(i=1:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  
  cp = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  cq = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  res = sapply(1:n_nodes, function(m) f_interval_approx_error(cp,cq,
                                                              boundaries[m,1],boundaries[m,2], 
                                                              Vectorize(demand.mixture.gamma),increasing = F,
                                                              l=1))
  list(t(res))
  
}
close(pb)
stopCluster(cl) 
toc()


save_path = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L1+dyadic_",n,"nodes",".RData",sep = "")
save(L1_Approx_Dyadic_50nodes,file=save_path)


i = sample(1:43848,1)
p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
p1 = p1[-1]

p2 = sort(unique(c(boundaries)))
p2 = p2[2:(length(p2)-1)]
q2 = L1_Approx_Dyadic_50nodes[[i]][,3]

sfun1 = stepfun(p1, q1, f = 0)
sfun2 = stepfun(p2, q2, f = 0)
plot(sfun1,xlim = c(0,1.1), ylim = c(0,1),
     ylab = "Quantity (MW)", xlab = paste('Price (', "\u20AC", '/MW)',sep=''),
     main = '',xaxs="i",yaxs="i", pch = 16)
lines(sfun2,col='blue',xlim = c(0,1.1),pch = 16)



#### 2.2 Marginal + L2 + dyadic ####

#### 2.2.1 Prior (n_nodes)/2 nodes placed by distribution and the remainder by dyadic ####
n_nodes = 50
n_prior_node = n_nodes/2 
probabilities= seq(0, 1, length.out = n_prior_node+1)[1:(n_prior_node)]
prior_node = quantile(day_ahead_demand$p_scaled,probabilities,type = 1)
prior_node = unname(prior_node)
prior_node

prior_l_boundaries = prior_node
prior_r_boundaries = c(prior_node[2:length(prior_node)],1)
prior_boundaries = cbind(prior_l_boundaries,prior_r_boundaries)

n_later_node = n_nodes/2


n_curve = max(day_ahead_demand$time)

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = n_curve
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

tic()
prior_error = foreach(cc=1:n_curve,.options.snow=opts)%dopar%{
  
  cp_real = day_ahead_demand[day_ahead_demand$time==cc,]$p_scaled
  cq_real = day_ahead_demand[day_ahead_demand$time==cc,]$q_cumsum_scaled
  Approx = sapply(1:n_prior_node, function(m) f_interval_approx_error(cp_real,cq_real,
                                                                      prior_boundaries[m,1],prior_boundaries[m,2], 
                                                                      Vectorize(demand.mixture.gamma),
                                                                      increasing = F,
                                                                      l=2))
  
}

close(pb)
stopCluster(cl)

Interval_E = cbind(Reduce("+", prior_error)[1,]/n_curve,
                   prior_l_boundaries,
                   prior_r_boundaries)
colnames(Interval_E) = c('error','l_boundary','r_boundary')
lst_CurvE = lapply(prior_error, function(x) {x = t(x)
                                             x = cbind(x,prior_l_boundaries,prior_r_boundaries)
                                             colnames(x) = c('error','approx_p','approx_q','l_boundary','r_boundary')
                                             x})
lst_Interval_E = list()


for(n in 1:n_later_node){
  
  p1 = unname(Interval_E[which.max(Interval_E[,1]),2])
  p3 = unname(Interval_E[which.max(Interval_E[,1]),3])
  p2 = unname((p1+p3)/2)
  
  Interval_drop = which.max(Interval_E[,1])
  Interval_E = Interval_E[-Interval_drop,]
    
  for(cc in 1:n_curve){

    print(paste0('n_later_node=',n,'; ','curve=',cc))
    E= lst_CurvE[[cc]]
    E = E[-Interval_drop,]
    
    cp = day_ahead_demand[day_ahead_demand$time==cc,]$p_scaled
    cq = day_ahead_demand[day_ahead_demand$time==cc,]$q_cumsum_scaled
    
    Il = c(p1,p2)
    Ir = c(p2,ifelse(p3 == max(cp),Inf,p3))
    
    
    left = f_interval_approx_error(cp,cq,Il[1],Il[2], Vectorize(demand.mixture.gamma),l=2,increasing = F)
    right = f_interval_approx_error(cp,cq,Ir[1],Ir[2], Vectorize(demand.mixture.gamma),l=2,increasing = F)
    
    
    E = rbind(E,
              c(left,p1,p2),
              c(right,p2,p3))
    
    lst_CurvE[[cc]] = E
    
  }
  
  
  Interval_mean_error = Reduce("+", lst_CurvE)[,1]/n_curve
  Interval_mean_error = tail(Interval_mean_error,2)
  Interval_E = rbind(Interval_E,
                     c(Interval_mean_error[1],p1,p2),
                     c(Interval_mean_error[2],p2,p3))
  lst_Interval_E[[n]] = Interval_E
  
}

toc()

save_path1 = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L2_lst_CurvE_",n_nodes,"nodes",".RData",sep = "")
save(lst_CurvE,file=save_path1)
save_path2 = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L2_lst_Interval_E_",n_nodes,"nodes",".RData",sep = "")
save(lst_Interval_E,file = save_path2)



##### 2.2.2. obtain the approximation with 50 nodes ####
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L2_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L2_lst_Interval_E_50nodes.RData")
n_node = 50
n = n_node/2
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
boundaries[,1] = real_boundaries[1:n_node]
boundaries[,2] = real_boundaries[2:(n_node+1)]

node = boundaries[,1]


cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

tic()
L2_Approx_Dyadic_50nodes = foreach(i=1:max(day_ahead_demand$time),.combine='rbind', .options.snow=opts,.packages = c("cubature","nor1mix")) %dopar% {
  p1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  q1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  
  numerator = sapply(1:(length(node)-1), function(n) f_CW(p1,q1,lower = node[-length(node)][n],upper = node[-1][n],f_weight=demand.mixture.gamma))
  denominator = sapply(1:(length(node)-1), function(n) f_W(lower = node[-length(node)][n],upper = node[-1][n],f_weight = demand.mixture.gamma))
  
  frac = numerator/denominator
  cn = f_CW(p1,q1,lower = node[length(node)-1], upper=node[length(node)], demand.mixture.gamma)/ f_W(lower = node[length(node)-1], upper=node[length(node)],demand.mixture.gamma)
  
  c_vector = c(frac, cn) - c(0,frac)           
  
  Approx = sapply(1:length(node), function(n) sapply(node, function(p) f_L2Approx(n,c_vector,node,p)))
  Approx = rowSums(Approx)
  Approx
}
close(pb)
stopCluster(cl) 
toc()

save_path = paste("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L2+dyadic_",n,"nodes2",".RData",sep = "")
save(L2_Approx_Dyadic_50nodes,file=save_path)

i = sample(1:43848,1)
cp1 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
cq1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled

cp2 = node
cq2 = L2_Approx_Dyadic_50nodes[i,]

cp1 = cp1[-1]
cp2 = cp2[-1]
sfun1 = stepfun(cp1, cq1, f = 0)
sfun2 = stepfun(cp2, cq2, f = 0)

plot(sfun1, xlim = c(0, max(cp1,cp2)), ylim = c(min(cq1, cq2), max(cq1, cq2)),
     ylab = "Scaled Quantity", xlab = "Scaled Price", 
     main = '',xaxs="i",yaxs="i",pch=16)
plot(sfun2,  col = "blue",  xlim = c(0, max(cp1,cp2)), add = TRUE,pch=16)

