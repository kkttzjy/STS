lambda2 = 0.1
c = 0.5
tau = 0.8
alpha_max = 1
maxit = 100
tree.num = 100
select.num = 100
SIS.num = 100 

X.low = t(eigenvectors(dm))[1:5,]

n = ncol(X.low)
p = nrow(X.low)

sigma = sd(X.low)
lambda1 = (sigma)^2*log(n)*p

### line search
#general numerical gradient function
grad.num = function(par, func, eps = 1e-6, ...) {
  par_new = par
  g = rep(0, length(par))
  f = func(par, ...)
  for (i in 1:length(par)) {
    par_new[i] = par[i] + eps
    f_new = func(par_new, ...)
    g[i] = (f_new - f) / eps
  }
  return (g)
}

#general line search function
line_search = function(par, func, grad = NULL, alpha_max = 1, c = 0.1, tau = 0.9, maxit = 100, diag = 0, ...) {
  f = func(par, ...)
  if (is.null(grad)) {
    g = grad.num(par, func, ...)
  } else {
    g = grad(par, ...)
  }
  if (diag > 0) {
    print("f")
    print(f)
    print("par")
    print(par)
    print("g")
    print(g)
  }
  alpha = alpha_max
  par_min = par
  f_min = f
  for (it in (1:maxit)) {
    par_new = par - alpha * g
    f_new = func(par_new, ...)
    if (diag > 1) print(c(f, f_new, alpha))
    if (f_new < f){
      f_min = f_new
      par_min = par_new
    }
    if (f - f_new > c * alpha * sum(g^2)) break
    alpha = tau * alpha
  }
  if (diag > 0) print(alpha)
  if (it == maxit) {
    if (f_min < f){
      return(par_min)
    }
    stop("line search failed")}
  return (par_new)
}


#tree score function for call back from line search
tree_score = function(par, X, Tree, lambda, i) {
  Tree_new = Tree
  Tree_new$Y[,i] = par
  proj = tree.project(X, Tree_new, lambda)
  f = sum(proj$D2)+sum(proj$D2.p)
  return (f)
}

#tree gradient function for call back from line search
tree_grad = function(par, X, Tree, lambda, i) {
  Tree_new = Tree
  Tree_new$Y[,i] = par
  proj = tree.project(X, Tree_new, lambda)
  g = gradient.calculation(X, Tree_new, proj, lambda, i)
  return (g)
}



################## Our algorithm functions #######################
#X: vector of length n
sum.squares = function(X) {
  return (sum(X^2))
}

N2E = function(N){
  N = t(N)
  N.max = max(N)
  E = matrix(F, N.max, N.max)
  E[N[, c(2, 1), drop = F]] = F
  E[N] = T
  E = t(E)
  E[N] = T
  return(E)
}

#X: matrix of p x n, coorindates of n data points
#Y1, Y2: vectors of length p, coordinates of two vertices
#return list(L, Phi, D)
#L: vector of length n, range in [0, 1], mapping of Phi on the edge from 0 to 1
#Phi: matrix of p x n, projection of X
#D2: vecotr of length n, squared distances
tree.edge.project = function(X, Y1, Y2, lambda2) {
  n = dim(X)[2]
  p = dim(X)[1]
  L = pmax(0, pmin(1, colSums((X - Y1) * (Y2 - Y1)) / sum.squares(Y2 - Y1)))
  Phi = Y1 %*% t(rep(1, n)) + (Y2 - Y1) %*% t(L)
  SS = 1/(n-1)*(X-Phi)%*%t(X-Phi)
  Sigma = matrix(0, p, p)
  diag(Sigma)= diag(SS)
  D2 = apply(X - Phi, 2, sum.squares)
  # X.mean = apply(X, 1, mean)
  # X.stand = apply(sweep(X,1,t(X.mean),"-"), 2,  sum.squares)
  # SD = sqrt(apply(sweep(X,1,t(X.mean),"-"), 2,  sum.squares)/(n-1))
  SD = sqrt(apply(X-Phi, 1, sum.squares)/n)
  # D2.SD =  1/2*apply((X-Phi)/SD, 2, sum.squares)
  D2.p = lambda2*sum.squares(Y1-Y2)
  return (list(L = L, Phi = Phi, D2 = D2, D2.p = D2.p, SD = SD))
}

#X: matrix of p x n, coorindates of n data points
#tree: Tree = list(Y, E, N)
#Tree$Y: matrix p x k, embedding of n vertices
#Tree$E: matrix of k x k, adjacency matrix, symmetric, diagonal all zeros
#Tree$N: matrtix of 2 x (k-1), vertice indices of edges
#return list(E, L, Phi, D2)
#I: matrix of 2 x n, indices of vertices for projected edges
#L: vector of length n, range in [0, 1], mapping of Phi on the edge from 0 to 1
#Phi: matrix of p x n, projection of X
#D2: vecotr of length n, squared distances
tree.project = function(X, Tree, lambda2) {
  p = dim(X)[1]
  n = dim(X)[2]
  k = dim(Tree$Y)[2]
  I = matrix(NA, 2, n)
  L = rep(NA, n)
  Phi = matrix(NA, p, n)
  D2 = rep(Inf, n)
  # D2.SD = rep(Inf, n)
  D2.p = rep(Inf, k-1)
  SD = matrix(NA, p, k-1)
  for (e in 1:(k-1)) {
    i = Tree$N[1, e]
    j = Tree$N[2, e]
    proj = tree.edge.project(X, Tree$Y[,i], Tree$Y[,j], lambda2)
    ind = proj$D2 < D2
    I[, ind] = c(i, j)
    L[ind] = proj$L[ind]
    Phi[, ind] = proj$Phi[, ind]
    D2[ind] = proj$D2[ind]
    # D2.SD[ind] = proj$D2.SD[ind]
    D2.p[e] = proj$D2.p
    SD[,e] = proj$SD
  }
  return(list(I = I, L = L, Phi = Phi, D2 = D2, D2.p = D2.p, SD = SD))
}

tree.update.node = function(X, Tree, proj, i, lambda2) {
  p = dim(X)[1]
  J = which(Tree$E[i,])
  m = length(J)
  numerator = 0
  denominator = 0
  for (j in J) {
    if (j > i) {
      t = (proj$I[1,] == i & proj$I[2,] == j)
      X1 = X[, t]
      L1 = 1 - proj$L[t]
      e = which(Tree$N[2,which(Tree$N[1,]==i)]==j)
    } else {
      t = (proj$I[1,] == j & proj$I[2,] == i)
      X1 = X[, t]
      L1 = proj$L[t]
      e = which(Tree$N[2,which(Tree$N[1,]==j)]==i)
    }
    
    
    sd = proj$SD[,e]
    if (length(L1) == 0){
      numpart = 0
    } else{
      if (length(L1) == 1){
        numpart = L1 * (X1 - Tree$Y[, j]*(1 - L1))
        
      } else{
        numpart = rowSums(sapply(seq(1,length(L1)), function(x){L1[x] * (X1[,x] - Tree$Y[, j]*(1 - L1[x]))}))
      }
    }
    numerator = numerator + numpart + lambda2*Tree$Y[, j]
    #numerator = numerator + rowSums(rep(L1, each = p) * (X1 - Tree$Y[, j] %*% t(1 - L1))) + lambda2*Tree$Y[, j]
    denominator = denominator + sum.squares(L1)
  }
  if (denominator == 0) {
    return(Tree$Y[,i])
  } else {
    return(numerator / (denominator+m*lambda2))
  }
}


#get edge list from adjancency matrix
E2N = function(E) {
  k = nrow(E)
  ind = which(E & upper.tri(E))
  N = rbind(((ind-1) %% k) + 1, floor((ind-1) / k) + 1)
  return(N)
}

initial.tree = function(p = 2, k = 2, X) {
  E = matrix(F, k, k)
  X.mean = apply(X, 1, mean)
  X.sd = apply(X, 1, sd)
  Y = matrix(rnorm(p*k, mean = X.mean, sd = X.sd), p, k)
  for (j in 1:k) {
    if (j > 1) {
      i = sample(1:(j-1), 1)
      E[i, j] = T
      E[j, i] = T
    }
  }
  perm = sample(1:k)
  E = E[perm, perm]
  return(list(Y = Y, E = E, N = E2N(E)))
}



SIS.initial = function(X, SIS.num, k = 2, maxit = 100, tol = 1e-2, lambda1 = 1, lambda2 = 1, select.num){
  p = nrow(X)
  n = ncol(X)
  X.stand = X
  tree.pool.current = NULL
  tree.pool.current.score = NULL
  tree.pool.current = lapply(1:SIS.num, function(x) initial.tree(p, k, X = X.stand))
  tree.pool.current.score = sapply(1:length(tree.pool.current), function(x){
    proj = tree.project(X.stand, tree.pool.current[[x]], lambda2)
    score = sum(proj$D2)+sum(proj$D2.p)
    score
  })
  
  for (m in 1: SIS.num){
    Tree = tree.pool.current[[m]]
    f_old = Inf
    highd.plot(Tree, X.stand)
    for (it in 1:maxit) {
      for (i in 1:k) {
        par = line_search(Tree$Y[,i], tree_score, NULL, c = c, tau = tau, X = X.stand, Tree = Tree, lambda = lambda2, i = i,diag =2)
        Tree$Y[,i] = par
        highd.plot(Tree, X.stand)
      }
      proj = tree.project(X.stand, Tree, lambda2)
      f = sum(proj$D2)+sum(proj$D2.p)
      if (f > f_old) stop("f > f_old")
      if (f_old - f < tol) break
      f_old = f
      proj = tree.project(X.stand, Tree, lambda2)
      Tree.distance = sapply(1:length(tree.pool.current), function(x)
        tree.disimilarity(Tree, tree.pool.current[[x]]))
      Tree.distance.min = min(Tree.distance)
      min.index = which(Tree.distance == Tree.distance.min)
      score = sum(proj$D2)+sum(proj$D2.p)
      if (score < tree.pool.current.score[min.index]){
        tree.pool.current[[min.index]] = Tree
        tree.pool.current.score[min.index] = score
      }
    }
    if (it == maxit) stop("maximum iteration reached")
      
      # if (Tree.distance.min > 5){
      #   score = sum(proj$D2)+sum(proj$D2.p)
      #   if (score < tree.pool.current.score[min.index]){
      #     tree.pool.current[[min.index]] = Tree
      #     tree.pool.current.score[min.index] = score
      #   }
      # }
    
  }
  score.seq.ord =  order(tree.pool.current.score, decreasing = F)
  select.tree = tree.pool.current[score.seq.ord[1:select.num]]
  return(select.tree)
}

tree.develop = function(old.tree, k, tree.num, X, lambda2, lambda3){
  n = ncol(X) 
  p = nrow(X) 
  X.mean = apply(X, 1, mean)
  X.sd = apply(X, 1, sd)
  new.tree = old.tree
  develop.event = sample(c(1,2), 1, prob = c(1-1/tree.num, 1/tree.num))
  if (develop.event == 1){
    v.select = sample(1:k, 1)
    index1 = old.tree$N[2,which(old.tree$N[1,]==v.select)]
    index2 = old.tree$N[1,which(old.tree$N[2,]==v.select)]
    index = c(index1, index2)
    ## calculate residuals
    proj = tree.project(X, old.tree, lambda2)
    residuals = proj$D2
    r.stand = residuals/(max(residuals))
    sampl.prob = exp(lambda3*r.stand)
    new.index = sample(seq(1,n,1), 1, prob = sampl.prob)
    y.new = X[,new.index]
    check.sim = apply(sweep(old.tree$Y, 1, y.new,"-"), 2, sum)
    
    if (any(check.sim == 0)){
      y.new = y.new + matrix(rnorm(p*1, mean = X.mean, sd = X.sd), p, 1)
    }
    
    new.tree$Y = cbind(old.tree$Y, y.new)
    new.tree$N = cbind(old.tree$N, c(v.select, k+1))
    new.tree$E = N2E(new.tree$N)
  } else{
    e.select.index = sample(seq(1,k-1,1), 1)
    e.select = old.tree$N[,e.select.index]
    y.new = (1/2)*(old.tree$Y[,e.select[1]]+old.tree$Y[,e.select[2]])
    new.tree$N = cbind(old.tree$N[,-e.select.index], c(e.select[1], k+1), c(e.select[2], k+1))
    new.tree$Y = cbind(old.tree$Y, y.new)
    new.tree$E = N2E(new.tree$N)
  }
  return(list(new.tree = new.tree, develop.event = develop.event))
}

SIS = function(X, SIS.num, tree.num, maxit = 100, tol = 1e-2, lambda1 = 1, lambda2 = 0.01, lambda3 = 1, select.num = 100, alpha = 1.5){
  #X.stand = t(scale(t(X)))
  X.stand = X
  #X.stand = t(pcareduceres)
  p = nrow(X.stand)
  n = ncol(X.stand)
  initial.tree.seq = SIS.initial(X.stand, SIS.num = SIS.num,  k = 2, select.num = select.num, lambda2 = lambda2)
  k=2
  old.tree.seq = initial.tree.seq
  Tree.diffsize.seq = initial.tree.seq[[1]]
  proj.initial = tree.project(X.stand, Tree.diffsize.seq, lambda2)
  score.diffsize.seq = sum(proj.initial$D2)+sum(proj.initial$D2.p)+lambda1*k^alpha
  score.new = score.min = score.diffsize.seq
  tree.pool.old = initial.tree.seq
  tree.current.index.his = NULL
  
  while( score.new <= score.min ){
    initial.index = seq(1:length(old.tree.seq))
    tree.current.index = initial.index
    tree.pool.current = lapply(1:length(old.tree.seq), function(x)
      tree.develop(old.tree.seq[[x]], k=k, tree.num = tree.num, X = X.stand, lambda2 = lambda2, lambda3 = lambda3)$new.tree)
    tree.pool.current.score1 = sapply(1:length(tree.pool.current), function(x){
      proj = tree.project(X.stand, tree.pool.current[[x]], lambda2)
      score = sum(proj$D2)+sum(proj$D2.p)
      score
    })
    # tree.pool.current.score2 = sapply(1:length(tree.pool.current), function(x){
    #   proj = tree.project(X.stand, tree.pool.current[[x]], lambda2)
    #   score = sum(proj$D2.SD)+sum(proj$D2.p)
    #   score
    # })
    for (m in 2: tree.num){
      new.tree.seq = lapply(1:length(old.tree.seq), function(x)
        tree.develop(old.tree.seq[[x]], k=k, tree.num = tree.num, X = X.stand, lambda2 = lambda2, lambda3 = lambda3)$new.tree)
      Tree.distance = sapply(1:length(new.tree.seq), function(y){
        score = sapply(1:length(tree.pool.current), function(x)
          tree.disimilarity(new.tree.seq[[y]], tree.pool.current[[x]]))
        score.min = min(score)[1]
        index.min = which(score == score.min)[1]
        res = t(c(score.min, index.min))
        return(res)
      })
      score.seq1 = sapply(1:length(new.tree.seq), function(x){
        proj = tree.project(X.stand, new.tree.seq[[x]], lambda2)
        score = sum(proj$D2)+sum(proj$D2.p)
        score
      })
      # score.seq2 = sapply(1:length(new.tree.seq), function(x){
      #   proj = tree.project(X.stand, new.tree.seq[[x]], lambda2)
      #   score = sum(proj$D2.SD)+sum(proj$D2.p)
      #   score
      # })
      score.pool.min = tree.pool.current.score1[t(Tree.distance[2,])]
      seq.current = t(rbind(Tree.distance, t(score.seq1), t(score.pool.min)))
      index1 = ifelse(seq.current[,1] > 1, 1, 0)
      index2 = ifelse(seq.current[,3]-seq.current[,4] < 0 , 1, 0)
      index3 = ifelse(index1 == 1 & index2 == 1, 1, 0)
      index4 = which(index3 == 1)
      tree.pool.current[seq.current[,2][index4]] = new.tree.seq[index4]
      tree.pool.current.score1[seq.current[,2][index4]] = score.seq1[index4]
      # tree.pool.current.score2[seq.current[,2][index4]] = score.seq2[index4]
      tree.current.index[seq.current[,2][index4]] = index4
    }
    tree.pool.current.initial = tree.pool.current
    tree.current.index.temp = tree.current.index
    
    for ( l in 1:length(tree.pool.current.initial)){
      Tree = tree.pool.current.initial[[l]]
      f_old = Inf
      for (it in 1:maxit) {
        for (i in 1:k) {
          par = line_search(Tree$Y[,i], tree_score, NULL, c = 0.1, tau = 0.8, X = X.stand, Tree = Tree, lambda = lambda2, i = i)
          Tree$Y[,i] = par
          #highd.plot(Tree, X.stand, trueplot = T, True.tree = True.tree)
        }
        proj = tree.project(X, Tree, lambda2)
        f = sum(proj$D2)+sum(proj$D2.p)
        if (f > f_old) stop("f > f_old")
        if (f_old - f < tol) break
        f_old = f
        proj = tree.project(X.stand, Tree, lambda2)
        Tree.distance = sapply(1:length(tree.pool.current), function(x)
          tree.disimilarity(Tree, tree.pool.current[[x]]))
        Tree.distance.min = min(Tree.distance)[1]
        min.index = which(Tree.distance == Tree.distance.min[1])[1]
        score1 = sum(proj$D2)+sum(proj$D2.p)
        # score2 = sum(proj$D2.SD)+sum(proj$D2.p)
        if (score1 < tree.pool.current.score1[min.index]){
          tree.pool.current[[min.index]] = Tree
          tree.pool.current.score1[min.index] = score1
          # tree.pool.current.score2[min.index] = score2
          tree.current.index[min.index] = tree.current.index.temp[l]
        }
        
        # if (Tree.distance.min > 0.4){
        #   score1 = sum(proj$D2)+sum(proj$D2.p)
        #   # score2 = sum(proj$D2.SD)+sum(proj$D2.p)
        #   if (score1 < tree.pool.current.score1[min.index]){
        #     tree.pool.current[[min.index]] = Tree
        #     tree.pool.current.score1[min.index] = score1
        #     # tree.pool.current.score2[min.index] = score2
        #     tree.current.index[min.index] = tree.current.index.temp[l]
        #   }
        # }
        
      }
      if (it == maxit) stop("maximum iteration reached")
      
    }

    score.seq.ord =  order(tree.pool.current.score1, decreasing = F)
    select.tree = tree.pool.current[score.seq.ord[1:select.num]]
    tree.current.index.ord = tree.current.index[score.seq.ord[1:select.num]]
    
    
    Tree.diffsize.seq[[k]] = select.tree[[1]]
    score.new = tree.pool.current.score1[score.seq.ord][1]+lambda1*(k+1)^alpha
    score.diffsize.seq = c(score.diffsize.seq, score.new)
    score.min = min(score.diffsize.seq)
    k = k + 1
    old.tree.seq = select.tree
    if (score.new <= score.min){
      tree.pool.old = append(tree.pool.old, select.tree)
      tree.current.index.his = cbind(tree.current.index.his, tree.pool.current.score1[tree.current.index.ord])
    }
  }
  opt.index = length(Tree.diffsize.seq)-1
  return(Tree = Tree.diffsize.seq[[opt.index]] )
}
###################### Tree simularity ########################
vertice.distance = function(V1, V2){
  n1 = length(V1[1,])
  n2 = length(V2[2,])
  p = length(V1[,1])
  distance = 0
  for ( i in 1:n1){
    distance_temp = sqrt(colSums((t(matrix(rep(V1[,i],each = n2),n2,p))-V2)^2))
    distance = distance + min(distance_temp)
  }
  return(distance)
}

edge.distance = function(T1, T2){
  N1 = T1$N
  N2 = T2$N
  n1 = dim(N1)[2]
  n2 = dim(N2)[2]
  V1 = T1$Y
  V2 = T2$Y
  p = length(V1[,1])
  D = 0
  for (i in 1:n1){
    index1 = N1[,i]
    V1.T1 = t(matrix(rep(V1[,index1[1]], each = n2), n2, p))
    V2.T1 = t(matrix(rep(V1[,index1[2]], each = n2), n2, p))
    V1.T2 = V2[,N2[1,],drop = F]
    V2.T2 = V2[,N2[2,],drop = F]
    D11 = sqrt(colSums((V1.T1-V1.T2)^2))
    D12 = sqrt(colSums((V1.T1-V2.T2)^2))
    D21 = sqrt(colSums((V2.T1-V1.T2)^2))
    D22 = sqrt(colSums((V2.T1-V2.T2)^2))
    D1 = apply(rbind(D11, D12), 2, min)
    D2 = apply(rbind(D21, D22), 2, min)
    D_temp = sum(sqrt(colSums(rbind(D1,D2)^2)))
    D = D + D_temp
  }
  return(D)
}

tree.disimilarity = function(T1, T2){
  disimilarity = vertice.distance(T1$Y, T2$Y)+edge.distance(T1,T2)
  return(disimilarity)
}

#################### Shorttest path algorithm ################################
### mat: mat[i,j] is the distance between nodes i and j, NA if i and j are not connected.

Edge.length = function(Y,N){
  V1 = N[1,]
  V2 = N[2,]
  length = sqrt(apply(Y[,V1]-Y[,V2], 2, sum.squares))
}

Adjacency.mat = function(Y,N){
  Elength = Edge.length(Y, N)
  k = dim(Y)[2]
  mat = matrix(NA,k,k)
  for ( i in 1:length(N[1,]) ){
    mat[N[1,i],N[2,i]] = Elength[i]
    mat[N[2,i],N[1,i]] = Elength[i]
  }
  return(mat)
}

#mat = Adjacency.mat(opt.tree1$Y, opt.tree1$N)
#ori = 1
spath = function (mat, ori) {
  k = dim(mat)[1]
  dist = rep(Inf, k)
  dist[ori] = 0
  updated = T
  while (updated) {
    updated = F
    for (i in 1:k) {
      for (j in 1:k) {
        if (!is.na(mat[i,j]) && dist[i] + mat[i,j] < dist[j]) {
          dist[j] = dist[i] + mat[i,j]
          updated = T
        }
      }
    }
  }
  return(dist)
}


highd.plot = function(Tree, X, proj = NULL, trueplot = F, True.tree = NULL){
  set.seed(12345)
  k = dim(Tree$Y)[2]
  n = dim(X)[2]
  PCA = prcomp(t(X), center = F, scale. = F)
  r = PCA$rotation
  #Tree.2d = t(t(Tree$Y)%*%r)
  Tree.2d = Tree$Y
  #xlimit = range(c(Tree.2d[1,], PCA$x[,1]))
  #ylimit = range(c(Tree.2d[2,], PCA$x[,2]))
  xlimit = range(c(Tree.2d[1,], X[1,]))
  ylimit = range(c(Tree.2d[2,], X[2,]))
  plot(t(Tree.2d[1:2,]), col = "blue",  pch = 19, lwd=3, xlim = xlimit, ylim = ylimit)
  for (e in 1:(k-1)) {
    lines(t(Tree.2d[, Tree$N[, e]]), lwd = 3, col = "black")
  }
  #points(PCA$x[,1],PCA$x[,2], col = "red", pch = 8)
  points(X[1,],X[2,], col = "red", pch = 8, lwd = 0.5)
  if (trueplot == T){
    True.tree.2d = t(t(True.tree$Y)%*%r)
    points(t(True.tree.2d[1:2,]), col = "orange", pch = 19)
    kt = dim(True.tree.2d)[2]
    for (e in 1:(kt-1)) {
      lines(t(True.tree.2d[, True.tree$N[, e]]), col = "green")
    }
  }
}

alpha = 1.5

proj.test1 = tree.project(X.stand, initial.tree.seq[[1]], lambda2)
score.size1= sum(proj.test1$D2)+sum(proj.test1$D2.p)+lambda1*2^alpha
proj.test2 = tree.project(X.stand, Tree.diffsize.seq[[2]], lambda2)
score.size2= sum(proj.test2$D2)+sum(proj.test2$D2.p)+lambda1*3^alpha
proj.test3= tree.project(X.stand, Tree.diffsize.seq[[3]], lambda2)
score.size3= sum(proj.test3$D2)+sum(proj.test3$D2.p)+lambda1*4^alpha
proj.test4= tree.project(X.stand, Tree.diffsize.seq[[4]], lambda2)
score.size4= sum(proj.test4$D2)+sum(proj.test4$D2.p)+lambda1*5^alpha
proj.test5= tree.project(X.stand, Tree.diffsize.seq[[5]], lambda2)
score.size5= sum(proj.test5$D2)+sum(proj.test5$D2.p)+lambda1*6^alpha


