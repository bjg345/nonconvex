source('functions.R')
load('finger_points.Rdata')

ind = 1:5000

start1  = Sys.time()
cov.select.fit1 = fit.water(grid.train[ind,], A.train[ind,ind], y=vals.train[ind], nu = .5, method='cov.select')
end1=Sys.time()
print(end1-start1)
start2  = Sys.time()
cov.select.fit2 = fit.water(grid.train[ind,], A.train[ind,ind], y=vals.train[ind], nu = .5, method='cov.select', n.neighbors = 15)
end2=Sys.time()
print(end2-start2)
save.image('out.Rdata')
