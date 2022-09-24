library(magrittr)
library(tidyverse)
z = qnorm(.975)
method.names = c("maxprec", "nearest", "weighted", "neighbor", "euclid", "BRISC")


out = matrix(NA, nrow=length(method.names), ncol = 2)
rownames(out) = method.names
colnames(out) = c('mse', 'coverage')

dat = data.frame(method = c(), n = c(), stdev = c(), mse = c(), coverage = c())

ind = c(1501:3000)[sapply(1501:3000, function(i) any(grepl(paste0('_', i, '.rds'), list.files('cs_fit') )) )]

getpar = function(j){

        if(j <= 1500){
                    n <- 250
                } else if (j <= 3000){
                    n <- 1200
                  } else (n = 10000)
        
                if((j %% 1500) <= 500){
                    stdev <- .1
                    } else if ((j %% 1500) <= 1000){
                        stdev <- .25
                    } else (stdev <- 1)
                    return(c(n, stdev))
            }

  for(method in method.names){
    
    
    err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2]; readRDS(file.path(method, paste0('err_', n, '_', stdev, '_', j, '.rds')))^2})
    mse = sapply(err2.list, mean) 
    
    covered = function(j){
       n=getpar(j)[1]; stdev=getpar(j)[2];   
      pred = readRDS(file.path(method, paste0('pred_', n, '_', stdev, '_', j, '.rds')))
      vals = readRDS(file.path(method, paste0('vals_', n, '_', stdev, '_', j, '.rds')))
      if(method != 'BRISC'){
        return(sapply(1:length(vals), function(i) vals[i] >= pred[[1,i]] - z*pred[[2,i]] & vals[i] <= pred[[1,i]] + z*pred[[2,i]]) %>% mean)
        }
      else{
        return(sapply(1:length(vals), function(i) vals[i] >= pred$prediction.ci[i,1] & vals[i] <= pred$prediction.ci[i,2]) %>% mean)
      }
    }
    
    coverage = sapply(ind, covered) 
    
    
    dat = rbind(dat, data.frame(method = method, n = sapply(ind, function(j) getpar(j)[1]), stdev=sapply(ind, function(j) getpar(j)[2]), mse = mse, coverage = coverage))
  }

out = dat %>% group_by(n, stdev, method) %>% summarise(mse = mean(mse), coverage = mean(coverage))

print(out)


out = matrix(NA, nrow=length(method.names), ncol = 2)
rownames(out) = method.names
colnames(out) = c('mse', 'coverage')

dat = data.frame(method = c(), n = c(), stdev = c(), mse = c(), coverage = c())

ind = c(1:1500)

getpar = function(j){

        if(j <= 1500){
                    n <- 250
                } else if (j <= 3000){
                    n <- 1200
                  } else (n = 10000)
        
                if((j %% 1500) <= 500){
                    stdev <- .1
                    } else if ((j %% 1500) <= 1000){
                        stdev <- .25
                    } else (stdev <- 1)
                    return(c(n, stdev))
            }

  for(method in method.names){
    
    
    err2.list = lapply(ind, function(j) {n=getpar(j)[1]; stdev=getpar(j)[2]; readRDS(file.path(method, paste0('err_', n, '_', stdev, '_', j, '.rds')))^2})
    mse = sapply(err2.list, mean) 
    
    covered = function(j){
       n=getpar(j)[1]; stdev=getpar(j)[2];   
      pred = readRDS(file.path(method, paste0('pred_', n, '_', stdev, '_', j, '.rds')))
      vals = readRDS(file.path(method, paste0('vals_', n, '_', stdev, '_', j, '.rds')))
      if(method != 'BRISC'){
        return(sapply(1:length(vals), function(i) vals[i] >= pred[[1,i]] - z*pred[[2,i]] & vals[i] <= pred[[1,i]] + z*pred[[2,i]]) %>% mean)
        }
      else{
        return(sapply(1:length(vals), function(i) vals[i] >= pred$prediction.ci[i,1] & vals[i] <= pred$prediction.ci[i,2]) %>% mean)
      }
    }
    
    coverage = sapply(ind, covered) 
    
    
    dat = rbind(dat, data.frame(method = method, n = sapply(ind, function(j) getpar(j)[1]), stdev=sapply(ind, function(j) getpar(j)[2]), mse = mse, coverage = coverage))
  }

out = dat %>% group_by(n, stdev, method) %>% summarise(mse = mean(mse), coverage = mean(coverage))

print(out)
