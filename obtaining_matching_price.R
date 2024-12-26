library(tidyverse)
library(foreach)
library(doParallel)
library(doSNOW)
library(pracma)

load("C:/Zehang_workspace/Approximation_Supply_Demand/match_prices.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/day_ahead_original_supply_curves.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/day_ahead_original_demand_curves.RData")

original_match_price = match_prices
day_ahead_supply = day_ahead_supply %>%
  group_by(date,hour) %>%
  mutate(time = cur_group_id()) %>%
  ungroup() %>%
  select(p = price,
         q = quantity,
         q_cumsum = cumsum_quantity,
         time)



day_ahead_demand = day_ahead_demand %>%
  group_by(date,hour) %>%
  mutate(time = cur_group_id()) %>%
  ungroup() %>%
  select(p = price,
         q = quantity,
         q_cumsum = cumsum_quantity,
         time)

max_q_cumsum = max(max(day_ahead_demand$q_cumsum),max(day_ahead_supply$q_cumsum))
max_p = 180.3

day_ahead_supply$p_scaled = day_ahead_supply$p/max_p
day_ahead_supply$q_cumsum_scaled = day_ahead_supply$q_cumsum / max_q_cumsum
day_ahead_supply$q_scaled =  day_ahead_supply$q / max_q_cumsum


day_ahead_demand$p_scaled = day_ahead_demand$p/max_p
day_ahead_demand$q_cumsum_scaled = day_ahead_demand$q_cumsum / max_q_cumsum
day_ahead_demand$q_scaled =  day_ahead_demand$q / max_q_cumsum 

##### 0. obtaining the matching prices from the original curves #####
cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

match_prices=foreach(i=1:43848,.combine='c', .options.snow=opts) %dopar% {  
  supply = day_ahead_supply[day_ahead_supply$time==i,]
  demand = day_ahead_demand[day_ahead_demand$time==i,]
  
  
  p1 = supply$p_scaled
  q1 = supply$q_cumsum_scaled
  p1 = p1[-1]
  
  p2 = demand$p_scaled
  q2 = demand$q_cumsum_scaled
  p2 = p2[-1]
  
  
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  
  difference1 = sfun1(p1) - sfun2(p1)
  WhereMatch1 =  which(difference1 >= 0)[1]
  match1 = p1[WhereMatch1]
  
  difference2 = sfun1(p2) - sfun2(p2)
  WhereMatch2 =  which(difference2 >= 0)[1]
  match2 = p2[WhereMatch2]
  
  pp = min(match1,match2)
  pp
}
close(pb)
stopCluster(cl)

match_prices 

#save(match_prices,file='C:/Zehang_workspace/Approximation_Supply_Demand/match_prices.RData')


i = 1
supply = day_ahead_supply %>%
  filter(time == i)
demand = day_ahead_demand %>%
  filter(time == i)

p1 = supply$p_scaled
q1 = supply$q_cumsum_scaled
p1 = p1[-1]

p2 = demand$p_scaled
q2 = demand$q_cumsum_scaled
p2 = p2[-1]


sfun1  <- stepfun(p1, q1, f = 0)
sfun2  <- stepfun(p2, q2, f = 0)


difference1 = sfun1(p1) - sfun2(p1)
WhereMatch1 =  which(difference1 >= 0)[1]
match1 = p1[WhereMatch1]

difference2 = sfun1(p2) - sfun2(p2)
WhereMatch2 =  which(difference2 >= 0)[1]
match2 = p2[WhereMatch2]

pp = min(match1,match2)


plot(sfun1, xlim = c(pp-0.01,pp +0.01),
     ylab = "Quantity", xlab = "Price",
     main = 'stepfun(p1, q1, f = 0)',xaxs="i",yaxs="i")
plot(sfun2,  col = "blue",  xlim = c(pp-0.01,pp +0.01), add = TRUE)
abline(v= pp,col='red')

##### 1. obtaining the matching prices from L2 + uniform approximation #####
n_node = 50
node= seq(0, 1, length.out = n_node)
supply_node = node
demand_node = node

load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_L2_Unif_GammaWeight_Approx_50nodes.RData")
supply_approx = Approx
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_L2_Unif_GammaWeight_Approx_50nodes.RData")
demand_approx = Approx

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

match_price=foreach(i=1:43848,.combine='c', .options.snow=opts) %dopar% {  
  supply = supply_approx[i,] 
  demand = demand_approx[i,]
  
  
  
  p1 = supply_node
  q1 = supply
  p1 = p1[-1]
  
  p2 = demand_node
  q2 = demand
  p2 = p2[-1]
  
  
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  
  difference1 = sfun1(p1) - sfun2(p1)
  WhereMatch1 =  which(difference1 >= 0)[1]
  match1 = p1[WhereMatch1]
  
  difference2 = sfun1(p2) - sfun2(p2)
  WhereMatch2 =  which(difference2 >= 0)[1]
  match2 = p2[WhereMatch2]
  
  pp = min(match1,match2)
  
  pp
}
close(pb)
stopCluster(cl)

mean(abs(original_match_price - match_price))*180.3
sqrt(mean((original_match_price - match_price)^2))*180.3

##### 2. obtaining the matching prices from L2 + marginal approximation #####
n_node = 50
probabilities= seq(0, 1, length.out = n_node)
supply_node = quantile(day_ahead_supply$p_scaled,probabilities,type = 1)
demand_node = quantile(day_ahead_demand$p_scaled,probabilities,type = 1)

load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/SP_Day-ahead_Supply_L2_Marginal_GammaWeight_Approx_50nodes.RData")
supply_approx = Approx
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_L2_marginal_GammaWeight_Approx_50nodes.RData")
demand_approx = Approx

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

match_price=foreach(i=1:43848,.combine='c', .options.snow=opts) %dopar% {  
  supply = supply_approx[i,] 
  demand = demand_approx[i,]
  
  
  
  p1 = supply_node
  q1 = supply
  p1 = p1[-1]
  
  p2 = demand_node
  q2 = demand
  p2 = p2[-1]
  
  
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  
  difference1 = sfun1(p1) - sfun2(p1)
  WhereMatch1 =  which(difference1 >= 0)[1]
  match1 = p1[WhereMatch1]
  
  difference2 = sfun1(p2) - sfun2(p2)
  WhereMatch2 =  which(difference2 >= 0)[1]
  match2 = p2[WhereMatch2]
  
  pp = min(match1,match2)
  
  pp
}
close(pb)
stopCluster(cl)

mean(abs(original_match_price - match_price))*180.3
sqrt(mean((original_match_price - match_price)^2))*180.3

##### 3. obtaining the matching prices from L2 + conditional approximation #####

n_node = 50
probabilities= seq(0, 1, length.out = n_node)
supply_Q =quantile(day_ahead_supply[(day_ahead_supply$p_scaled != 0),]$q_scaled,0.125)
demand_Q =quantile(day_ahead_demand[(day_ahead_demand$p_scaled != 0),]$q_scaled,0.125)
supply_node = quantile(day_ahead_supply[day_ahead_supply$q_scaled > supply_Q,]$p_scaled,probabilities,type = 1)
demand_node = quantile(day_ahead_demand[day_ahead_demand$q_scaled > demand_Q,]$p_scaled,probabilities,type = 1)

load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/SP_Day-ahead_Supply_L2_Conditional_GammaWeight_Approx_50nodes.RData")
supply_approx = Approx
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_L2_conditional_GammaWeight_Approx_50nodes.RData")
demand_approx = Approx

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

match_price=foreach(i=1:43848,.combine='c', .options.snow=opts) %dopar% {  
  supply = supply_approx[i,] 
  demand = demand_approx[i,]
  
  
  
  p1 = supply_node
  q1 = supply
  p1 = p1[-1]
  
  p2 = demand_node
  q2 = demand
  p2 = p2[-1]
  
  
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  
  difference1 = sfun1(p1) - sfun2(p1)
  WhereMatch1 =  which(difference1 >= 0)[1]
  match1 = p1[WhereMatch1]
  
  difference2 = sfun1(p2) - sfun2(p2)
  WhereMatch2 =  which(difference2 >= 0)[1]
  match2 = p2[WhereMatch2]
  
  pp = min(match1,match2)
  
  pp
}
close(pb)
stopCluster(cl)

mean(abs(original_match_price - match_price))*180.3
sqrt(mean((original_match_price - match_price)^2))*180.3


##### 4. obtaining the matching prices from L2 + dyadic approximation #####

load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_L2+dyadic_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_L2_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_L2_lst_Interval_E_50nodes.RData")
n = 50

supply_approx = L2_Approx_Dyadic_50nodes

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

supply_node = boundaries[,1]

load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_L2+dyadic_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_L2_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_L2_lst_Interval_E_50nodes.RData")
n = 50

demand_approx = L2_Approx_Dyadic_50nodes

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

demand_node = boundaries[,1]


cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

# curves at 9385th hour don't cross
match_price=foreach(i=1:43848,.combine='c', .options.snow=opts) %dopar% {  
  supply = supply_approx[i,] 
  demand = demand_approx[i,] 
  
  
  
  p1 = supply_node
  q1 = supply
  p1 = p1[-1]
  
  p2 = demand_node
  q2 = demand
  p2 = p2[-1]
  
  
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  
  difference1 = sfun1(p1) - sfun2(p1)
  WhereMatch1 =  which(difference1 >= 0)[1]
  match1 = p1[WhereMatch1]
  
  difference2 = sfun1(p2) - sfun2(p2)
  WhereMatch2 =  which(difference2 >= 0)[1]
  match2 = p2[WhereMatch2]
  
  pp = min(match1,match2)
  
  pp =ifelse(is.na(pp),max(supply_node,demand_node),pp)
}
close(pb)
stopCluster(cl)

mean(abs(original_match_price - match_price))*180.3
sqrt(mean((original_match_price - match_price)^2))*180.3



##### 5. obtaining the matching prices from L2 + Mixed nodes approximation #####
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L2_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L2_lst_Interval_E_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L2+dyadic_25nodes.RData")


supply_approx = L2_Approx_Dyadic_50nodes

n_nodes = 50
n = n_nodes/2
n_curve = max(day_ahead_supply$time)
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
boundaries[,1] = real_boundaries[1:n_nodes]
boundaries[,2] = real_boundaries[2:(n_nodes+1)]

supply_node = boundaries[,1]

load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L2_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L2_lst_Interval_E_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L2+dyadic_25nodes.RData")


demand_approx = L2_Approx_Dyadic_50nodes

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


demand_node = boundaries[,1]


cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

# curves at 9385th hour don't cross
match_price=foreach(i=1:43848,.combine='c', .options.snow=opts) %dopar% {  
  supply = supply_approx[i,] 
  demand = demand_approx[i,] 
  
  
  
  p1 = supply_node
  q1 = supply
  p1 = p1[-1]
  
  p2 = demand_node
  q2 = demand
  p2 = p2[-1]
  
  
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  
  difference1 = sfun1(p1) - sfun2(p1)
  WhereMatch1 =  which(difference1 >= 0)[1]
  match1 = p1[WhereMatch1]
  
  difference2 = sfun1(p2) - sfun2(p2)
  WhereMatch2 =  which(difference2 >= 0)[1]
  match2 = p2[WhereMatch2]
  
  pp = min(match1,match2)
  
  #pp =ifelse(is.na(pp),max(supply_node,demand_node),pp)
}
close(pb)
stopCluster(cl)

mean(abs(original_match_price - match_price))*180.3
sqrt(mean((original_match_price - match_price)^2))*180.3



##### 6. obtaining the matching prices from L1 + uniform approximation #####

n_node = 50
node= seq(0, 1, length.out = n_node+1)
boundaries = cbind(node[1:50],node[2:51])
node = boundaries[,1]
supply_node = node
demand_node = node

load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_L1_Unif_GammaWeight_Approx_50nodes.RData")
supply_approx = L1_Approx_Unif_50nodes
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_L1_Unif_GammaWeight_Approx_50nodes.RData")
demand_approx = L1_Approx_Unif_50nodes

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

match_price=foreach(i=1:43848,.combine='c', .options.snow=opts) %dopar% {  
  supply = supply_approx[[i]][,3]
  demand = demand_approx[[i]][,3]
  
  
  
  p1 = supply_node
  q1 = supply
  p1 = p1[-1]
  
  p2 = demand_node
  q2 = demand
  p2 = p2[-1]
  
  
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  
  difference1 = sfun1(p1) - sfun2(p1)
  WhereMatch1 =  which(difference1 >= 0)[1]
  match1 = p1[WhereMatch1]
  
  difference2 = sfun1(p2) - sfun2(p2)
  WhereMatch2 =  which(difference2 >= 0)[1]
  match2 = p2[WhereMatch2]
  
  pp = min(match1,match2)
  
  pp
}
close(pb)
stopCluster(cl)

mean(abs(original_match_price - match_price))*180.3
sqrt(mean((original_match_price - match_price)^2))*180.3


##### 7. obtaining the matching prices from L1 + margianl approximation #####

n_node = 50
probabilities= seq(0, 1, length.out = n_node+1)
node = quantile(day_ahead_supply$p_scaled,probabilities,type = 1)
boundaries = cbind(node[1:50],node[2:51])
supply_node = boundaries[,1]

node = quantile(day_ahead_demand$p_scaled,probabilities,type = 1)
boundaries = cbind(node[1:50],node[2:51])
demand_node = boundaries[,1]


load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_L1_marginal_GammaWeight_Approx_50nodes.RData")
supply_approx = L1_Approx_marginal_50nodes
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_L1_marginal_GammaWeight_Approx_50nodes.RData")
demand_approx = L1_Approx_marginal_50nodes

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

match_price=foreach(i=1:43848,.combine='c', .options.snow=opts) %dopar% {  
  supply = supply_approx[[i]][,3]
  demand = demand_approx[[i]][,3]
  
  
  
  p1 = supply_node
  q1 = supply
  p1 = p1[-1]
  
  p2 = demand_node
  q2 = demand
  p2 = p2[-1]
  
  
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  
  difference1 = sfun1(p1) - sfun2(p1)
  WhereMatch1 =  which(difference1 >= 0)[1]
  match1 = p1[WhereMatch1]
  
  difference2 = sfun1(p2) - sfun2(p2)
  WhereMatch2 =  which(difference2 >= 0)[1]
  match2 = p2[WhereMatch2]
  
  pp = min(match1,match2)
  
  pp
}
close(pb)
stopCluster(cl)

mean(abs(original_match_price - match_price))*180.3
sqrt(mean((original_match_price - match_price)^2))*180.3


##### 8. obtaining the matching prices from L1 + conditional approximation #####

Q =quantile(day_ahead_supply[(day_ahead_supply$p_scaled != 0),]$q_scaled,0.125)
# Q is the 12.5% quantile of the quantitis that are included in non-zero offers

Q = unname(Q);Q
n_node = 50
probabilities= seq(0, 1, length.out = n_node+1)
node = quantile(day_ahead_supply[day_ahead_supply$q_scaled > Q,]$p_scaled,probabilities,type = 1)
boundaries = cbind(node[1:50],node[2:51])
supply_node = boundaries[,1]


Q =quantile(day_ahead_demand[(day_ahead_demand$p_scaled != 0),]$q_scaled,0.125)
# Q is the 12.5% quantile of the quantitis that are included in non-zero offers

Q = unname(Q);Q

n_node = 50
probabilities= seq(0, 1, length.out = n_node+1)
node = quantile(day_ahead_demand[day_ahead_demand$q_scaled > Q,]$p_scaled,probabilities,type = 1)
boundaries = cbind(node[1:50],node[2:51])
demand_node = boundaries[,1]


load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_L1_conditional_GammaWeight_Approx_50nodes.RData")
supply_approx = L1_Approx_conditional_50nodes
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_L1_conditional_GammaWeight_Approx_50nodes.RData")
demand_approx = L1_Approx_conditional_50nodes

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

match_price=foreach(i=1:43848,.combine='c', .options.snow=opts) %dopar% {  
  supply = supply_approx[[i]][,3]
  demand = demand_approx[[i]][,3]
  
  
  
  p1 = supply_node
  q1 = supply
  p1 = p1[-1]
  
  p2 = demand_node
  q2 = demand
  p2 = p2[-1]
  
  
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  
  difference1 = sfun1(p1) - sfun2(p1)
  WhereMatch1 =  which(difference1 >= 0)[1]
  match1 = p1[WhereMatch1]
  
  difference2 = sfun1(p2) - sfun2(p2)
  WhereMatch2 =  which(difference2 >= 0)[1]
  match2 = p2[WhereMatch2]
  
  pp = min(match1,match2)
  
  pp
}
close(pb)
stopCluster(cl)

mean(abs(original_match_price - match_price))*180.3
sqrt(mean((original_match_price - match_price)^2))*180.3


##### 9. obtaining the matching prices from L1 + dyadic approximation #####
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_L1+dyadic_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_L1_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_L1_lst_Interval_E_50nodes.RData")
n = 50

supply_approx = L1_Approx_Dyadic_50nodes

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

supply_node = boundaries[,1]

load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_L1_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_L1_lst_Interval_E_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_L1+dyadic_50nodes.RData")
n = 50

demand_approx = L1_Approx_Dyadic_50nodes

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

demand_node = boundaries[,1]

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

match_price=foreach(i=1:43848,.combine='c', .options.snow=opts) %dopar% {  
  supply = supply_approx[[i]][,3]
  demand = demand_approx[[i]][,3]
  
  
  
  p1 = supply_node
  q1 = supply
  p1 = p1[-1]
  
  p2 = demand_node
  q2 = demand
  p2 = p2[-1]
  
  
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  
  difference1 = sfun1(p1) - sfun2(p1)
  WhereMatch1 =  which(difference1 >= 0)[1]
  match1 = p1[WhereMatch1]
  
  difference2 = sfun1(p2) - sfun2(p2)
  WhereMatch2 =  which(difference2 >= 0)[1]
  match2 = p2[WhereMatch2]
  
  pp = min(match1,match2)
  
  pp
}
close(pb)
stopCluster(cl)

mean(abs(original_match_price - match_price))*180.3
sqrt(mean((original_match_price - match_price)^2))*180.3


##### 10. obtaining the matching prices from L1 + Mixed Nodes approximation #####


load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L1_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L1_lst_Interval_E_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_GammaWeight_MixedNode_L1+dyadic_25nodes.RData")


supply_approx = L1_Approx_Dyadic_50nodes

n_nodes = 50
n = n_nodes/2
n_curve = max(day_ahead_supply$time)
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
boundaries[,1] = real_boundaries[1:n_nodes]
boundaries[,2] = real_boundaries[2:(n_nodes+1)]

supply_node = boundaries[,1]


load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L1_lst_CurvE_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L1_lst_Interval_E_50nodes.RData")
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_GammaWeight_MixedNode_L1+dyadic_25nodes.RData")


demand_approx = L1_Approx_Dyadic_50nodes

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

demand_node = boundaries[,1]

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

match_price=foreach(i=1:43848,.combine='c', .options.snow=opts) %dopar% {  
  supply = supply_approx[[i]][,3]
  demand = demand_approx[[i]][,3]
  
  
  
  p1 = supply_node
  q1 = supply
  p1 = p1[-1]
  
  p2 = demand_node
  q2 = demand
  p2 = p2[-1]
  
  
  sfun1  <- stepfun(p1, q1, f = 0)
  sfun2  <- stepfun(p2, q2, f = 0)
  
  
  difference1 = sfun1(p1) - sfun2(p1)
  WhereMatch1 =  which(difference1 >= 0)[1]
  match1 = p1[WhereMatch1]
  
  difference2 = sfun1(p2) - sfun2(p2)
  WhereMatch2 =  which(difference2 >= 0)[1]
  match2 = p2[WhereMatch2]
  
  pp = min(match1,match2)
  
  pp
}
close(pb)
stopCluster(cl)

mean(abs(original_match_price - match_price))*180.3
sqrt(mean((original_match_price - match_price)^2))*180.3



##### 11. obtaining the matching prices from RBF approximation #####
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_Approx_mesh_free_M_50.RData")
l_matrix_coeff_supply = Approx
load("C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_Approx_mesh_free_M_50.RData")
l_matrix_coeff_demand = Approx

M_demand = 50
M_supply = 50



cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = 43848
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
#
match_price=foreach(i=1:43848,.combine='c', .options.snow=opts,.packages = c('pracma')) %dopar% {
  matrix_coeff_supply = l_matrix_coeff_supply[[i]]
  matrix_coeff_demand = l_matrix_coeff_demand[[i]]
  
  
  supply_cq1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  f_supply_iteration_at_nodes = function(x){sapply(1:M_supply, function(n) {matrix_coeff_supply[n,1]*( erf((x-matrix_coeff_supply[n,2] )/ matrix_coeff_supply[n,3])+1)})}
  f_supply_approx = function(x) {res=apply(f_supply_iteration_at_nodes(x),MARGIN = 1,sum);return(res+supply_cq1[1])}#;return(res+cq1[length(cq1)])
  
  
  demand_cq1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  demand_cq1 = rev(demand_cq1)
  matrix_coeff_demand = Approx[[i]]
  f_demand_iteration_at_nodes = function(x){sapply(1:M_demand, function(n) {matrix_coeff_demand[n,1]*( erf((x-matrix_coeff_demand[n,2] )/ matrix_coeff_demand[n,3])+1)})}
  f_demand_approx = function(x) {res=apply(f_demand_iteration_at_nodes(x),MARGIN = 1,sum);return(res+demand_cq1[1])}
  x = seq(1e-8,1,length.out = 1e5)

  difference = f_supply_approx(x) - rev(f_demand_approx(x))
  idx = which.min(abs(difference))
  x[idx]
  #match_price =c(match_price,x[idx])
}
close(pb)
stopCluster(cl)

mean(abs(original_match_price - match_price))*180.3
sqrt(mean((original_match_price - match_price)^2))*180.3

i = 2
matrix_coeff_supply = l_matrix_coeff_supply[[i]]
matrix_coeff_demand = l_matrix_coeff_demand[[i]]

###########
supply_cq1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
f_supply_iteration_at_nodes = function(x){sapply(1:M_supply, function(n) {matrix_coeff_supply[n,1]*( erf((x-matrix_coeff_supply[n,2] )/ matrix_coeff_supply[n,3])+1)})}
f_supply_approx = function(x) {res=apply(f_supply_iteration_at_nodes(x),MARGIN = 1,sum);return(res+supply_cq1[1])}#;return(res+cq1[length(cq1)])


demand_cq1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
demand_cq1 = rev(demand_cq1)
f_demand_iteration_at_nodes = function(x){sapply(1:M_demand, function(n) {matrix_coeff_demand[n,1]*( erf((x-matrix_coeff_demand[n,2] )/ matrix_coeff_demand[n,3])+1)})}
f_demand_approx = function(x) {res=apply(f_demand_iteration_at_nodes(x),MARGIN = 1,sum);return(res+demand_cq1[1])}

x = seq(0,1,by = 0.001)
plot(x,f_supply_approx(x),ylim=c(0,1))
lines(x,rev(f_demand_approx(x)),ylim=c(0,1))

cp3 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
cq3 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
cp3 = cp3[-1]
sfun3 = stepfun(cp3, cq3, f = 0)

cp4 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
cq4 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
cp4 = cp4[-1]
sfun4 = stepfun(cp4, cq4, f = 0)

lines(sfun3,ylim=c(0,1))
lines(sfun4,ylim=c(0,1))
