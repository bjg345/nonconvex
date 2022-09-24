library(BRISC)

if(!dir.exists('weighted')) dir.create('weighted')
id = as.numeric(commandArgs(trailingOnly = T))
if(id <= 1500){
  n = 250
} else if (id <= 3000){
  n = 1200
} else (n = 10000)

if((id %% 1500) <= 500){
  stdev = .1
} else if ((id %% 1500) <= 1000){
  stdev = .25
} else (stdev = 1)

#if(file.exists(file.path('weighted', paste0('vals_', n, '_', stdev, '_', id, '.rds')))) quit()

load('../finger_data.Rdata')
source('../functions.R')

id = as.numeric(commandArgs(trailingOnly = T))
if(id <= 1500){
  n = 250
} else if (id <= 3000){
  n = 1200
} else (n = 10000)

if((id %% 1500) <= 500){
  stdev = .1
} else if ((id %% 1500) <= 1000){
  stdev = .25
} else (stdev = 1)
set.seed(id)

ind.train = sample.int(10000, .8*n, replace = F)
noise.train = rnorm(.8*n, sd = stdev)

ind.test = sample.int(10000, .2*n, replace = F)
noise.test = rnorm(.2*n, sd = stdev)

grid.train = grid.train[ind.train,]
grid.test = grid.test[ind.test,]

ground.train = ground.train[ind.train]
ground.test = ground.test[ind.test]

vals.train = ground.train + noise.train
vals.test = ground.test + noise.test

A.train = A.train[ind.train, ind.train]
A.test = A.test[ind.test, ind.test]

if(n > 250) k = 10 else k = NULL
if(n <= 250) {
	weighted.fit = fit.water(grid.train, A.train, y=vals.train, nu = .5, method='cov.select', n.neighbors=k)
    }else{
		weighted.fit = readRDS(file.path('cs_fit', paste0('fit_', n, '_', stdev, '_', id, '.rds')))
	} 
weighted.pred = sapply(1:nrow(grid.test), function(i) pred_neighbor(new.loc = grid.test[i,], A=A.train, A.vec = NULL, 
                                          grid.train, vals.train, weighted.fit$par, nu=.5, method = 'precweighted', n.neighbors = 10, 
                                           D=as.matrix(dist(grid.train)), water) )

weighted.err = unlist(weighted.pred[1,]) - ground.test

saveRDS(weighted.fit, file.path('weighted', paste0('fit_', n, '_', stdev, '_', id, '.rds')))
saveRDS(weighted.pred, file.path('weighted', paste0('pred_', n, '_', stdev, '_', id, '.rds')))
saveRDS(weighted.err, file.path('weighted', paste0('err_', n, '_', stdev, '_', id, '.rds')))
saveRDS(vals.test, file.path('weighted', paste0('vals_', n, '_', stdev, '_', id, '.rds')))

