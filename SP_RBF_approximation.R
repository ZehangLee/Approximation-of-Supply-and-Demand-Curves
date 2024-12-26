library(vroom)
library(tidyverse)
library(tidyselect)
library(pracma)
library(doSNOW)
library(foreach)
library(doParallel)
library(tictoc)

library(mixtools)
library(LaplacesDemon)
library(nor1mix)
library(cubature)

load("C:/Zehang_workspace/Approximation_Supply_Demand/gammamix_match_prices.RData")

########################################################################
# onestepdata
########################################################################

onestepdata = function (q_low_bound,q_high_bound,h,p,q_cumsum){
  q1 = q_cumsum
  p1 = p
  i = 1
  
  while(q1[i] <= q_low_bound){
    q1[i] = q_low_bound
    i = i +1
  }
  I1 = i - 1 # locate the prices lower than the prices where the curve jumps
  
  while( (q1[i] < q_high_bound) & (i < length(q1)) ){
    i = i +1
  }
  I2 = i -1  # locate the prices higher than the prices where the curve jumps
  
  while(i <= length(q1)){
    q1[i] = q_high_bound
    i = i +1
  }
  
  
  m = round(p1[length(p1)]/h)
  datax = seq(0, m*h, by= h)
  dataf = matrix(0, length(datax), 1)
  i = 1
  for(n in 1:length(p1)){
    while ( (i <= m) & (datax[i] <= p1[n]) ){
      dataf[i] = q1[n]
      i = i +1
    }
  }
  
  dataf[i] = q_high_bound
  center = round( (p1[I1]/h + (p1[I2] - p[I2])/h) * .5 )
  
  return(list('datax'=datax,
              'dataf'=dataf,
              'center'=center))
}







########################################################################
# RBF
########################################################################

RBF_approx_coeff = function(p,q,M_supply,shape = 1/3000){
  q_cumsum = cumsum(q)
  
  an_matrix = matrix(0, nrow= 4,ncol = M_supply +1)
  # the first row is prices at M steps
  # the second row is quantities at M steps
  # the third row is offers can buy when reaching the M prices
  # the fourth row is number of steps before reaching the M prices
  an_matrix[1, M_supply+1] = p[length(p)]
  an_matrix[2, M_supply+1] = q[length(q)]
  an_matrix[3, M_supply+1] = q_cumsum[length(q_cumsum)]
  an_matrix[4, M_supply+1] = 1 # initialize to be 1, it means need a step to reach the final prices
  an_matrix[4, 1] = 1 # initialize to be 1, it means need a step to jump from 0
  
  iter1 = 1
  for(num in 0:(M_supply-1)){
    qcumsum_jump_scale = an_matrix[3,M_supply-num+1]/(M_supply-num )
    j = 0
    
    # while( (iter1 <= length(q_cumsum))  &  
    #        (q_cumsum[length(q_cumsum) - iter1] > (an_matrix[3, M_supply -num +1] -qcumsum_jump_scale) )){
    #   iter1 = iter1 + 1
    #   j = j + 1
    # }
    # 
    # an_matrix[1, M_supply - num] = p[length(p) - iter1]
    # an_matrix[2, M_supply - num] = q[length(q) - iter1]
    # an_matrix[3, M_supply - num] = q_cumsum[length(q_cumsum) - iter1]
    # an_matrix[4, M_supply - num +1] = an_matrix[4, M_supply - num +1] + j
    
    tryCatch(
      {
        while( (iter1 <= length(q_cumsum))  &
               (q_cumsum[length(q_cumsum) - iter1] > (an_matrix[3, M_supply -num +1] -qcumsum_jump_scale) )){
          iter1 = iter1 + 1
          j = j + 1
        }
        
        an_matrix[1, M_supply - num] = p[length(p) - iter1]
        an_matrix[2, M_supply - num] = q[length(q) - iter1]
        an_matrix[3, M_supply - num] = q_cumsum[length(q_cumsum) - iter1]
        an_matrix[4, M_supply - num +1] = an_matrix[4, M_supply - num +1] + j
      },
      
      error=function(error_message){})
    
  }
  
  
  matrix_coeff = matrix(0,M_supply, 3)
  
  for(i in 1:M_supply){
    q_cumsum1 = an_matrix[3,i]
    q_cumsum2 = an_matrix[3,i+1]
    jump_num = an_matrix[4, i+1]
    center = an_matrix[1, i]
    
    if(jump_num == 1){
      matrix_coeff[i,] = c((q_cumsum2 - q_cumsum1)/2, center, 0.5)
    }
    
    if(jump_num > 1){
      h = 0.001
      paras = onestepdata(q_cumsum1, q_cumsum2, h, p, q_cumsum)
      datax = paras$datax
      dataf = paras$dataf
      center = paras$center
      a1 = (q_cumsum2 -q_cumsum1)/2
      a2 = center
      a3 = shape #2000 # 3000
      a4 = q_cumsum1
      z0 = c(a2, a3)
      
      F = function(z, zdata) a1 * ( erf( z[2] * (datax- h*z[1]) ) +1) + a4
      z = lsqcurvefit(F,z0,datax,dataf)$x #optimization
      matrix_coeff[i,] = c(a1,h*z[1],1/z[2])
    }
  }
  return(matrix_coeff)
}


#### weight function being used in the SUPPLY CURVE approximation procedure ####
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


#### weight function being used in the DEMAND CURVE approximation procedure ####
dgammamixEM <- function(x,GammaDist) {
  dist = sapply(1:length(GammaDist$lambda),function(i) GammaDist$lambda[i]*dgamma(x, GammaDist$gamma.pars[1,i], 
                                                                                  1/GammaDist$gamma.pars[2,i]))
  return(sum(dist))
}

demand.mixture.gamma =  function(x) dgammamixEM(x,gammamixEM_p_m2)






############################################
# Approximation of supply curve of Spanish day-ahead market
###########################################

load("C:/Zehang_workspace/Approximation_Supply_Demand/day_ahead_original_supply_curves.RData")

max_q_cumsum = 90834.4#max(max(day_ahead_demand$q_cumsum),max(day_ahead_supply$q_cumsum))
max_p = 180.3


day_ahead_supply = day_ahead_supply %>%
  group_by(date,hour) %>%
  mutate(time = cur_group_id()) %>%
  ungroup() %>%
  select(p = price,
         q = quantity,
         q_cumsum = cumsum_quantity,
         time)
day_ahead_supply$p_scaled = day_ahead_supply$p/max_p
day_ahead_supply$q_cumsum_scaled = day_ahead_supply$q_cumsum / max_q_cumsum
day_ahead_supply$q_scaled =  day_ahead_supply$q / max_q_cumsum



M_supply = 50


cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
#
tic()
Approx = foreach(i=1:max(day_ahead_supply$time), .options.snow=opts,.packages = c("pracma")) %dopar% {
  p = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  q = day_ahead_supply[day_ahead_supply$time==i,]$q_scaled
  
  matrix_coeff_supply = RBF_approx_coeff(p,q,M_supply)
}
close(pb)
stopCluster(cl) 
toc()

save(Approx,file='C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_Approx_mesh_free_M_50.RData')

i = sample(1:16,1)
cp1 = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
cq1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled


matrix_coeff_supply = Approx[[i]]
f_iteration_at_nodes = function(x){sapply(1:M_supply, function(n) {matrix_coeff_supply[n,1]*( erf((x-matrix_coeff_supply[n,2] )/ matrix_coeff_supply[n,3])+1)})}
f_approx = function(x) {res=apply(f_iteration_at_nodes(x),MARGIN = 1,sum);return(res+cq1[1])}#;return(res+cq1[length(cq1)])
x= seq(0,1,by =0.01)
plot(x,f_approx(x))

cp1 = cp1[-1]
sfun1 = stepfun(cp1, cq1, f = 0)

plot(sfun1, ylim=c(0,1),
     ylab = "Scaled Quantity", xlab = "Scaled Price", 
     main = '',xaxs="i",yaxs="i",pch=16)
lines(x,f_approx(x),  col = "blue", ylim=c(0,1), pch=16)





load('C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_supply_Approx_mesh_free_M_50.RData')
cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

#
EmprDis_err=foreach(i=25:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("pracma","nor1mix")) %dopar% {
  p = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  p1 = p[-1]
  sfun1  = stepfun(p1, q1, f = 0)
  
  matrix_coeff_supply = Approx[[i]]
  f_iteration_at_nodes = function(x){sapply(1:M_supply, function(n) {matrix_coeff_supply[n,1]*( erf((x-matrix_coeff_supply[n,2] )/ matrix_coeff_supply[n,3])+1)})}
  f_approx = function(x) {res=apply(f_iteration_at_nodes(x),MARGIN = 1,sum);return(res+q1[1])} 
  
  
  f = function(x) (f_approx(x) - sfun1(x))^2 *supply.mixture.gamma(x)
  
  dist <- integrate(f,
                    upper = Inf,
                    lower = 0,
                    #subdivisions=3000, 
                    rel.tol=.Machine$double.eps^.05)
  dist$value
  
  
}
close(pb)
stopCluster(cl)

mean_EmprDis_err_nodes=mean(EmprDis_err)
sqrt(mean_EmprDis_err_nodes)

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_supply$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
EmprDis_err=foreach(i=25:max(day_ahead_supply$time),.combine='c', .options.snow=opts,.packages = c("pracma","nor1mix")) %dopar% {
  p = day_ahead_supply[day_ahead_supply$time==i,]$p_scaled
  q1 = day_ahead_supply[day_ahead_supply$time==i,]$q_cumsum_scaled
  p1 = p[-1]
  sfun1  = stepfun(p1, q1, f = 0)
  
  matrix_coeff_supply = Approx[[i]]
  f_iteration_at_nodes = function(x){sapply(1:M_supply, function(n) {matrix_coeff_supply[n,1]*( erf((x-matrix_coeff_supply[n,2] )/ matrix_coeff_supply[n,3])+1)})}
  f_approx = function(x) {res=apply(f_iteration_at_nodes(x),MARGIN = 1,sum);return(res+q1[1])} 
  
  
  
  f = function(x) abs(f_approx(x) - sfun1(x)) *supply.mixture.gamma(x)
  
  dist <- integrate(f,
                    upper = Inf,
                    lower = 0,
                    #subdivisions=3000, 
                    rel.tol=.Machine$double.eps^.05)
  dist$value
  
  
}
close(pb)
stopCluster(cl)

mean_EmprDis_err_nodes=mean(EmprDis_err)
mean_EmprDis_err_nodes




############################################
# Approximation of demand curve of Spanish day-ahead market
###########################################

load("C:/Zehang_workspace/Approximation_Supply_Demand/day_ahead_original_demand_curves.RData")

max_q_cumsum = 90834.4#max(max(day_ahead_demand$q_cumsum),max(day_ahead_supply$q_cumsum))
max_p = 180.3


day_ahead_demand = day_ahead_demand %>%
  group_by(date,hour) %>%
  mutate(time = cur_group_id()) %>%
  ungroup() %>%
  select(p = price,
         q = quantity,
         q_cumsum = cumsum_quantity,
         time)

day_ahead_demand$p_scaled = day_ahead_demand$p/max_p
day_ahead_demand$q_cumsum_scaled = day_ahead_demand$q_cumsum / max_q_cumsum
day_ahead_demand$q_scaled =  day_ahead_demand$q / max_q_cumsum


M_demand = 50


cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
#
tic()
Approx = foreach(i=1:max(day_ahead_demand$time), .options.snow=opts,.packages = c("pracma")) %dopar% {
  p = 1-day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  p = rev(p)
  q = day_ahead_demand[day_ahead_demand$time==i,]$q_scaled
  q = rev(q)
  
  
  matrix_coeff_demand = RBF_approx_coeff(p,q,M_demand,shape = 1/500)
}
close(pb)
stopCluster(cl) 
toc()

save(Approx,file='C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_Approx_mesh_free_M_50.RData')

i = sample(1:43848,1)
cp1 = 1 - day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
cp1 = rev(cp1)
cq1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
cq1 = rev(cq1)


matrix_coeff_demand = Approx[[i]]
f_iteration_at_nodes = function(x){sapply(1:M_demand, function(n) {matrix_coeff_demand[n,1]*( erf((x-matrix_coeff_demand[n,2] )/ matrix_coeff_demand[n,3])+1)})}
f_approx = function(x) {res=apply(f_iteration_at_nodes(x),MARGIN = 1,sum);return(res+cq1[1])}#;return(res+cq1[length(cq1)])

cp2 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
cq2 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
cp2 = cp2[-1]
sfun2 = stepfun(cp2, cq2, f = 0)
x = seq(0,1,by=0.001)
plot(sfun2, ylim=c(0,1),
     ylab = "Scaled Quantity", xlab = "Scaled Price", 
     main = '',xaxs="i",yaxs="i",pch=16)
lines(x,rev(f_approx(x)),  col = "red", ylim=c(0,1), pch=16)




load('C:/Zehang_workspace/Approximation_Supply_Demand/gamma weight approximation/sp_day-ahead_demand_Approx_mesh_free_M_50.RData')
cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

#
EmprDis_err=foreach(i=25:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("pracma","nor1mix")) %dopar% {
  cp1 = 1 - day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  cp1 = rev(cp1)
  cq1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  cq1 = rev(cq1)

  
  matrix_coeff_demand = Approx[[i]]
  f_iteration_at_nodes = function(x){sapply(1:M_demand, function(n) {matrix_coeff_demand[n,1]*( erf((x-matrix_coeff_demand[n,2] )/ matrix_coeff_demand[n,3])+1)})}
  f_approx = function(x) {res=apply(f_iteration_at_nodes(x),MARGIN = 1,sum);return(res+cq1[1])}#;return(res+cq1[length(cq1)])
  
  cp2 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  cq2 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  cp2 = cp2[-1]
  sfun2 = stepfun(cp2, cq2, f = 0)
  
  f = function(x) (rev(f_approx(x)) - sfun2(x))^2 *demand.mixture.gamma(x)
  
  dist <- integrate(f,
                    upper = Inf,
                    lower = 0,
                    #subdivisions=3000, 
                    rel.tol=.Machine$double.eps^.05)
  dist$value
  
  
}
close(pb)
stopCluster(cl)

mean_EmprDis_err_nodes=mean(EmprDis_err)
rmse = sqrt(mean_EmprDis_err_nodes)

cl = makeSOCKcluster(16)
registerDoSNOW(cl)
iterations = max(day_ahead_demand$time) -25 +1
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
EmprDis_err=foreach(i=25:max(day_ahead_demand$time),.combine='c', .options.snow=opts,.packages = c("pracma","nor1mix")) %dopar% {
  cp1 = 1 - day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  cp1 = rev(cp1)
  cq1 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  cq1 = rev(cq1)

  
  matrix_coeff_demand = Approx[[i]]
  f_iteration_at_nodes = function(x){sapply(1:M_demand, function(n) {matrix_coeff_demand[n,1]*( erf((x-matrix_coeff_demand[n,2] )/ matrix_coeff_demand[n,3])+1)})}
  f_approx = function(x) {res=apply(f_iteration_at_nodes(x),MARGIN = 1,sum);return(res+cq1[1])}#;return(res+cq1[length(cq1)])
  
  cp2 = day_ahead_demand[day_ahead_demand$time==i,]$p_scaled
  cq2 = day_ahead_demand[day_ahead_demand$time==i,]$q_cumsum_scaled
  cp2 = cp2[-1]
  sfun2 = stepfun(cp2, cq2, f = 0)
  
  f = function(x) abs(rev(f_approx(x)) - sfun2(x)) *demand.mixture.gamma(x)
  
  dist <- integrate(f,
                    upper = Inf,
                    lower = 0,
                    #subdivisions=3000, 
                    rel.tol=.Machine$double.eps^.05)
  dist$value
  
  
}
close(pb)
stopCluster(cl)

mean_EmprDis_err_nodes=mean(EmprDis_err)
mae=mean_EmprDis_err_nodes


