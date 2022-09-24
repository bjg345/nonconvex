neighbor.loglik = function(param, nu, y, dists, coords, neighbor.list, ord){
  
  range = exp(param[1])
  nugget = exp(param[2])
  sigma = exp(param[3])
  
  n = length(y)
  
  ind.loglik = function(i){
    
    neighbors = neighbor.list[[i]]
    
    if(is.null(neighbor.list[[i]])) {return(dnorm(y[i], mean=0, sd = sqrt(sigma^2+nugget), log=T))}
    
    
    
    
    
    ind = c(i, neighbors)
    
    
    dist.nn = dists[ind, ind]
    if(any(dist.nn == Inf)){
      #if there are neighbors not connected in water, start from nearest neighbor and move outward without adding unconnected points
      ind = i
      
      neighbor.ord = sort(dists[i, neighbors], index.return=T)$ix
      for(j in neighbor.ord){
        ind.prop = c(ind, neighbors[j])
        if(!any(dists[neighbors[j], ind.prop]==Inf)) ind = ind.prop
      }
      
      neighbors = setdiff(ind, i)
      dist.nn = dists[ind, ind]
    }
    
    cov.mat =geostatsp::matern(dist.nn, param = c(range=range, variance=sigma^2, shape = nu)) + 
      diag(nugget, length(ind)) 
    
    Sig12 = cov.mat[1, 2:length(ind), drop=F]
    Sig22inv = solve(cov.mat[2:length(ind), 2:length(ind)])
    
    mu.cond = Sig12 %*% Sig22inv %*% y[neighbors]
    var.cond = cov.mat[1,1] - Sig12 %*% Sig22inv %*% t(Sig12)
    
    return(dnorm(y[i], mean = mu.cond, sd=sqrt(var.cond), log=T) )
    
  }
  
  return(-sum(sapply(1:n, ind.loglik)))
  
}


cov.select.loglik = function(param, nu = nu, y, dists, rip, n.neighbors = NULL, coords){
  
  range = exp(param[1])
  nugget = exp(param[2])
  sigma = exp(param[3])
  print(c(range,nugget,sigma))
  n = length(y)
  
  cliques = sapply(rip$cliques, as.numeric)
  seps = sapply(rip$separators, as.numeric)
  
  if(!is.null(n.neighbors)){
    require(BRISC2)
    
    clique.lik = function(clique){
      
      if(length(clique)==0) return(0)
    
      coords = coords[clique,]
      y = y[clique]
      
      res = BRISC2_estimation(coords = coords, y = y, cov.model = 'exponential', n.neighbors = min(n.neighbors, length(clique)),
      phi = 1/range, tau.sq = nugget, sigma.sq = sigma^2, verbose=F)$log_likelihood
      
      return(res)
    }
    
    out = -Reduce('+', lapply(cliques, clique.lik)) + Reduce('+', lapply(seps, clique.lik))
    
    if(is.nan(out) | out == Inf | out == -Inf) out = 9e99
    print(out)
    return(out)
    
  }
  
  else{  
    out = -Reduce('+', lapply(cliques, function(clique){
      cov.mat = geostatsp::matern(dists[clique, clique], param = c(range=range, variance=sigma^2, shape = nu)) + 
        diag(nugget, length(clique))
      
      return(Rfast::dmvnorm(y[clique], mu = rep(0, length(clique)), sigma = cov.mat, logged=T))
      
    })) + Reduce('+', lapply(seps, function(sep){
      
      cov.mat = geostatsp::matern(dists[sep, sep], param = c(range=range, variance=sigma^2, shape = nu)) + 
        diag(nugget, length(sep))
      
      return(Rfast::dmvnorm(y[sep], mu = rep(0, length(sep)), sigma = cov.mat, logged=T))
      
    }))
    
    if(is.nan(out) | out == Inf | out == -Inf) out = 1e50
    
    return(out)
  }
}


pred_neighbor = function(new.loc, grid, y, param, nu = .5, n.neighbors, A, A.vec = NULL, D, water, thresh = Inf, method = NULL){
  require(igraph)
  require(fields)
  
  range = exp(param[1])
  nugget = exp(param[2])
  sigma = exp(param[3])
  
  euc.dists = fields::rdist(matrix(new.loc, nrow=1), grid)
  
  
  if(is.null(A.vec)){
    require(rgeos)
    within.thresh = which(euc.dists < thresh)
    
    l <- vector("list", length(within.thresh))
    
    for (j in seq_along(within.thresh)) {
      l[[j]] <- Lines(list(Line(as.matrix(rbind(new.loc, grid[within.thresh[j],])))), as.character(within.thresh[j]))
    }
    
    con = rep(F, nrow(grid))
    con[within.thresh] = gCovers(water, SpatialLines(l), byid=T)
  }
  else{con = A.vec}
  
  adj.dists = euc.dists/con
  
  cutoff = sort(adj.dists)[n.neighbors]
  if(cutoff == Inf) {
    neighbor.ind = which(adj.dists < Inf)
  } else{
    neighbor.ind = which(adj.dists <= cutoff)[1:n.neighbors]
  }
  
  if(length(neighbor.ind) == 0) return(list(mu=NA, sd = NA))
  
  if(method == 'euclidean'){
    
    dist.nn = D[neighbor.ind, neighbor.ind]
    dist.nn = rbind(adj.dists[neighbor.ind], dist.nn)
    dist.nn = cbind(c(0, adj.dists[neighbor.ind]), dist.nn)
    
    L = nrow(dist.nn)
    
    cov.mat =geostatsp::matern(dist.nn, param = c(range=range, variance=sigma^2, shape = nu)) + 
      diag(nugget, L)
    
    Sig11 = cov.mat[1,1]
    Sig12 = cov.mat[1, 2:L, drop=F]
    Sig22inv = solve(cov.mat[2:L, 2:L])
    
    mu.post = Sig12 %*% Sig22inv %*% y[neighbor.ind]
    sd.post =sqrt( Sig11 - Sig12 %*% Sig22inv %*% t(Sig12) )
    
    return(list(mu = mu.post, sd = sd.post))
    
  }
  
  neighbor.adj = A[neighbor.ind, neighbor.ind]
  L = length(neighbor.ind)
  
  if(method == 'nearest.clique'){
    
    if(!all(neighbor.adj + diag(L) == 1)){
      
      neighbor.ind.red = c()
      
      ord = neighbor.ind[sort(adj.dists[neighbor.ind], index.return=T)$ix]
      for(candidate in ord){
        
        if(all(A[c(neighbor.ind.red, candidate), c(neighbor.ind.red, candidate)] +
               diag(1+length(neighbor.ind.red) )== 1)){
          neighbor.ind.red = c(neighbor.ind.red, candidate)
        }
        
      }
      
    } else{
      neighbor.ind.red = neighbor.ind
    }
    
    dist.nn = D[neighbor.ind.red, neighbor.ind.red]
    dist.nn = rbind(adj.dists[neighbor.ind.red], dist.nn)
    dist.nn = cbind(c(0, adj.dists[neighbor.ind.red]), dist.nn)
    
    L = nrow(dist.nn)
    
    cov.mat =geostatsp::matern(dist.nn, param = c(range=range, variance=sigma^2, shape = nu)) + 
      diag(nugget, L)
    
    Sig11 = cov.mat[1,1]
    Sig12 = cov.mat[1, 2:L, drop=F]
    Sig22inv = solve(cov.mat[2:L, 2:L])
    
    mu.post = Sig12 %*% Sig22inv %*% y[neighbor.ind.red]
    sd.post =sqrt( Sig11 - Sig12 %*% Sig22inv %*% t(Sig12) )
    
    return(list(mu = mu.post, sd = sd.post))
  }
  
  
  g = graph_from_adjacency_matrix(neighbor.adj, 'undirected')
  g = set.vertex.attribute(g, "name", value=1:length(neighbor.ind))
  
  if(method == 'maxprec') {cliques.list = max_cliques(g)}
  else if(method == 'precweighted'){
    i = 1
    cliques.list = list()
    while(vcount(g)>0){
      clique = attr(largest_cliques(g)[[1]], 'name')
      cliques.list[[i]] = as.numeric(clique)
      g = delete.vertices(g, clique)
      
      i = i+1
    }
  }  
  
  
  
  out = lapply(cliques.list, function(clique){
    
    neighbor.ind.clique = neighbor.ind[clique]
    
    dist.nn = D[neighbor.ind.clique, neighbor.ind.clique]
    dist.nn = rbind(adj.dists[neighbor.ind.clique], dist.nn)
    dist.nn = cbind(c(0, adj.dists[neighbor.ind.clique]), dist.nn)
    
    L = nrow(dist.nn)
    
    cov.mat =geostatsp::matern(dist.nn, param = c(range=range, variance=sigma^2, shape = nu)) + 
      diag(nugget, L)
    
    Sig11 = cov.mat[1,1]
    Sig12 = cov.mat[1, 2:L, drop=F]
    Sig22inv = solve(cov.mat[2:L, 2:L])
    
    mu.post = Sig12 %*% Sig22inv %*% y[neighbor.ind.clique]
    prec.post = 1/(Sig11 - Sig12%*%Sig22inv%*%t(Sig12))
    
    return(list(mu.post=mu.post, prec.post=prec.post))
  }
  )
  
  mu.posts = sapply(out, function(x) x$mu.post)
  prec.posts = sapply(out, function(x) x$prec.post)
  
  if (method == 'precweighted') return(list(mu=sum(mu.posts*prec.posts)/sum(prec.posts), sd = (sum(prec.posts))^(-.5)))
  else if(method == 'maxprec') return(list(mu=mu.posts[which.max(prec.posts)], sd = (max(prec.posts))^(-.5)))
}


fit.water = function(coords, adj, y, nu = .5, method, n.neighbors=NULL){
  
  if(!method %in% c('neighbor', 'cov.select')) stop('invalid method name')
  
  if(method == 'cov.select'){
    require(igraph)
    require(gRbase)
    require(BRISC)
    
    D = as.matrix(stats::dist(coords))
    
    g.unchord = graph_from_adjacency_matrix(adj, 'undirected')
    g.chord = is_chordal(g.unchord, newgraph=T)$newgraph
    rip = rip(as(g.chord, 'graphNEL'))
    
    BRISC.fit = BRISC_estimation(y=y, coords=coords)
    range.start = 1/BRISC.fit$Theta[3]
    nugget.start = BRISC.fit$Theta[2]
    sigma.start = sqrt(BRISC.fit$Theta[1])
    
    fit = optim(par = log(c(range.start,nugget.start,sigma.start)), fn = cov.select.loglik, nu = nu, y = y, dists = D, rip = rip,
                 n.neighbors = n.neighbors, coords = coords)
    
    return(fit)
  } else if(method =='neighbor'){
    require(igraph)
    
    D = as.matrix(stats::dist(coords))
    D.adj = D/adj
    diag(D.adj) = 0
    
    g = graph_from_adjacency_matrix(adj, 'undirected') 
    ord = sort(rowSums(coords), index.return=T)$ix
    
    neighbor.list = lapply(1:nrow(coords), function(i) {
      
      if(i==ord[[1]]){ return(NULL)}
      
      dists = D.adj[i,]
      
      eligible = ord[1:(-1+(which(ord==i)))]
      #points eligible to be neighbors need to be previous in graph ordering
      
      eligible.dists = dists[eligible]
      if(all (eligible.dists==Inf)) return(NULL)
      
      if(n.neighbors > length(eligible)){
        return(eligible[which(eligible.dists != Inf)]) 
      }  
      
      else{
        cutoff = sort(eligible.dists)[n.neighbors] #find what distance kth nearest neighbor is
        if(cutoff == Inf) return(eligible[which(eligible.dists<Inf)])
        return(eligible[which(eligible.dists <= cutoff)][1:n.neighbors]) #return neighbors within that distance (final part is in case of ties)
      }
      
    }
    )
    
    fit = optim(par = c(.5,.5,.5), fn = neighbor.loglik, nu = nu, neighbor.list = neighbor.list, y = y, dists = D.adj, coords = coords,
                 ord = ord)
    
    return(list(fit = fit, neighbor.list = neighbor.list))
    
  }
  
}
