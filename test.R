ind = sample(25000, 1000)

load('u_data.Rdata')

A = A[ind,ind]

grid = grid[ind,]

ground = ground[ind]

save(list=c('A', 'grid', 'water', 'ground'), file='u_samp.Rdata')


ind = sample(25000, 1000)

load('finger_data.Rdata')

A = A[ind,ind]

grid = grid[ind,]

ground = ground[ind]

save(list=c('A', 'grid', 'water', 'ground'), file='finger_samp.Rdata')




