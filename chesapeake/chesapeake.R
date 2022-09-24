library(tidyverse)
library(magrittr)
library(sf)
library(blockCV)
library(viridis)
library(BRISC)
library(MLmetrics)
library(rgeos)
library(ggplot2)
library(data.table)
library(raster)
library(gdistance)

nu = .5

chesa = st_read("mygeodata/chesapeake-polygon.shp")
chesa.detail = st_read("Chesapeake_Bay_Shoreline_Medium_Resolution/Chesapeake_Bay_Shoreline_Medium_Resolution.shp")
projcrs <- st_crs(chesa)

#x = extent(chesa)
#buff.rec = 10000
#x@xmin=x@xmin-buff.rec; x@xmax = x@xmax+buff.rec; x@ymin = x@ymin-buff.rec; x@ymax = x@ymax+buff.rec
#ras <- raster(nrow = 1000, ncol = 1000, ext = x)
#values(ras) <- 1
#ras <- mask(ras, chesa)
#tr <- transition(ras, transitionFunction=mean, directions=8) %>% geoCorrection('c')

#data is all TWQM dissolved oxygon 1/1/2010 - 12/31/2019
ph.dat= fread('WaterQualityWaterQualityStation.csv', quote="") %>%
  filter(Parameter =="\"PH\""  ) %>%
  group_by(Station) %>%
  summarise(y = mean(MeasureValue, na.rm=T),
            lat = mean(Latitude), lon = mean(Longitude)) %>% #checked no variation in lat/lon
  dplyr::select(lat, lon, y) 

ph.sf <- st_as_sf(ph.dat, coords = c('lon', 'lat'), crs = 'WGS84')

water <- as_Spatial(chesa)# %>% spTransform(CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
coords <- SpatialPoints(ph.dat[, c('lon', 'lat')],
                        proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')) %>% 
  spTransform(projcrs$wkt)




#dists <- costDistance(tr, coords@coords, coords@coords)
#diag(dists) = NA


con = function(i, p.grid) {
  require(sp)
  
  n = nrow(p.grid)
  
  l <- vector("list", i)
  
  for (j in 1:i) {
    l[[j]] <- Lines(list(Line(as.matrix(rbind(p.grid[i,], p.grid[j,])))), as.character(j))
  }
  
  out = gCovers(water, 
                SpatialLines(l, proj4string=coords@proj4string), byid=T)
  print(i)
  return(out)
  
}

out.list = sapply(1:nrow(ph.dat), FUN = con, p.grid = coords@coords)
A = matrix(nrow = nrow(ph.dat), ncol=nrow(ph.dat))
for(i in 1:nrow(A)){
  A[i, 1:i] = out.list[[i]]
}
A[upper.tri(A)] = t(A)[upper.tri(A)]
diag(A) = 0
hist(rowSums(A), breaks = 15)


saveRDS(ph.dat, 'ph.rds')

p1=ggplot() + 
  ggtitle("Chesapeake pH") +
  geom_sf(data = chesa.detail, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = y), data = ph.sf) +
  coord_sf() +
  scale_color_viridis(name='pH')
ggsave('pH.png', p1)

p2=ggplot() + 
  ggtitle("Chesapeake pH - buffered boundary") +
  geom_sf(data = chesa, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = y), data = ph.sf) +
  coord_sf() +
  scale_color_viridis(name='pH')
ggsave('PH_buff.png', p2)



n <- nrow(ph.dat)

set.seed(89)
blocked =spatialBlock(coords, rows=10, cols=10, selection='checkerboard') #blocked test/train desgin
blocked$foldID[69]=2
blocked$foldID[136]=1

train.id = which(blocked$foldID==1)
test.id = which(blocked$foldID==2)

train.dat <- ph.dat[train.id,]
test.dat <- ph.dat[test.id,]

#for(i in 1:nrow(A)){
  
 # neigh <- base::sort(dists[i,train.id], na.last = T, index.return=T)$ix[1:3]
 # A[i, train.id[neigh]] <- A[train.id[neigh], i] <- 1
  
  
#}



coords.norm <- scale(coords@coords)
BRISC.mod = BRISC_estimation(coords = coords.norm[train.id,], y = train.dat$y, cov.model = 'exponential',
                             n.neighbors=15)
BRISC.pred = BRISC_prediction(BRISC.mod, coords.0 = coords.norm[test.id,])




cs.fit = fit.water(coords.norm[train.id,], A[train.id, train.id], y= ph.dat$y[train.id] ,method = 'cov.select')
neighbor.fit = fit.water(coords.norm[train.id,], A[train.id, train.id], y= ph.dat$y[train.id] ,method = 'neighbor',
                         n.neighbors=15)

cs.pred = sapply(test.id, function(i)
  pred_neighbor(new.loc = coords.norm[i,], A=A[train.id,train.id], A.vec = A[i, train.id], 
                coords.norm[train.id,], ph.dat$y[train.id], cs.fit$par, nu=.5, 
                method = 'nearest.clique', n.neighbors = 15,
                D=as.matrix(dist(coords.norm[train.id, ]))) )
cs.pred.vals = cs.pred[1,] %>% unlist

neighbor.pred = sapply(test.id, function(i)
  pred_neighbor(new.loc = coords.norm[i,], A=A[train.id,train.id], A.vec = A[i, train.id], 
                coords.norm[train.id,], ph.dat$y[train.id], neighbor.fit$fit$par, nu=.5, 
                method = 'euclidean', n.neighbors = 15,
                D=as.matrix(dist(coords.norm[train.id, ]))) )
neighbor.pred.vals = neighbor.pred[1,] %>% unlist

BRISC.err <- BRISC.pred$prediction - test.dat$y
BRISC.err.sf <- st_as_sf(ph.dat[test.id,] %>% mutate(err = BRISC.err)
                         , coords = c('lon', 'lat'), crs = 'WGS84')

p3=ggplot() + 
  ggtitle("BRISC error") +
  geom_sf(data = chesa, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = err), data = BRISC.err.sf) +
  coord_sf() +
  scale_color_viridis(limits=c(-.3, 1.11))
ggsave('BRISC_error.png', p3)

cs.err <- cs.pred.vals - test.dat$y
cs.err.sf <- st_as_sf(ph.dat[test.id,] %>% mutate(err = cs.err)
                         , coords = c('lon', 'lat'), crs = 'WGS84')
p4=ggplot() + 
  ggtitle("Covariance selection error") +
  geom_sf(data = chesa, size = .3, color = "black", fill = "cyan1") + 
  geom_sf(mapping = aes(color = err), data = cs.err.sf) +
  coord_sf() +
  scale_color_viridis(limits=c(-.3,1.11))
ggsave('cs_error.png', p4)

MSE(BRISC.pred$prediction, ph.dat$y[test.id])
MSE(cs.pred.vals, ph.dat$y[test.id])
MSE(neighbor.pred.vals, ph.dat$y[test.id])
MSE(mean(ph.dat$y[test.id]), ph.dat$y[test.id])
