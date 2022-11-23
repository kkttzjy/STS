#################### Shorttest path algorithm ################################
## linear tree
Edge.length = function(Y,N){
  V1 = N[1,]
  V2 = N[2,]
  length = sqrt(apply(Y[,V1]-Y[,V2], 2, sum.squares))
}

Adjacency.mat.linear = function(Y,N){
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
spath.linear = function (mat, ori) {
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


## curved tree
### mat: mat[i,j] is the distance between nodes i and j, NA if i and j are not connected.
arc.length.partial = function(cell, V, alpha){
    knots = seq(0,1,length.out = 21)  
    nknots = length(knots)
    f.knots = matrix(unlist(lapply(seq(1,nknots),function(k) pc(V,cell,alpha,knots[k]))), nrow=p, ncol=nknots)
    f1 = f.knots[,-1]
    f2 = f.knots[,-21]
    f.piece = sqrt(apply(f1-f2, 2, sum.squares))
    length = sum(f.piece)
  return(length)
}

arc.length = function(G){
  V1 = G$N[1,]
  V2 = G$N[2,]
  length = rep(0, length(G$N[1,]))
  for ( i in 1:length(G$N[1,]) ){
  Y1 = G$Y[,V1[i]]
  Y2 = G$Y[,V2[i]]
  knots = seq(0,1,length.out = 21)  
  nknots = length(knots)
  f.knots = matrix(unlist(lapply(seq(1,nknots),function(k) pc(Y1,Y2,G$alpha[,i],knots[k]))), nrow=p, ncol=nknots)
  f1 = f.knots[,-1]
  f2 = f.knots[,-21]
  f.piece = sqrt(apply(f1-f2, 2, sum.squares))
  length[i] = sum(f.piece)
  }
  return(length)
}
# arc.length = function(b,Y1,Y2,alpha){
#   integrand = function(x) {sqrt(sum.squares(Y1+x*(Y2-Y1)+x*(x-1)*alpha))}
#   value = integrate(integrand, lower = 0, upper = b)$value
#   
#   return(s)
# }

# Edge.length = function(Y,N){
#   V1 = N[1,]
#   V2 = N[2,]
#   length = sqrt(apply(Y[,V1]-Y[,V2], 2, sum.squares))
# }

Adjacency.mat = function(G){
  Elength = arc.length(G)
  k = dim(G$Y)[2]
  mat = matrix(NA,k,k)
  for ( i in 1:length(G$N[1,]) ){
    mat[G$N[1,i],G$N[2,i]] = Elength[i]
    mat[G$N[2,i],G$N[1,i]] = Elength[i]
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

### pseudotime calculation
## pseudotime assignment
opt.G = G
proj.opt = tree.curved.project(X, opt.G, lambda2)

### order data according to the true tree
true.proj  = tree.curved.project(X, True.G, lambda2)
## pseudotime assignment
# initial cell
inital.cell.true.index = sample(seq(1:k),1)

mat0 = Adjacency.mat(True.G)
node.dis0 = spath(mat0, inital.cell.true.index)

# pseudotime calculation based on projection
n = dim(X)[2]
index.edge.true = true.proj$I
cell.pt.true = rep(0,n)
for ( i in 1:n){
  e = true.proj$E[i]
  E = index.edge.true[,i]
  cell = true.proj$fi[,i]
  V1 = True.G$Y[,E[1]]
  V2 = True.G$Y[,E[2]]
  alpha = True.G$alpha[,e]
  if ( node.dis0[E[1]] < node.dis0[E[2]]){
    V.old = V1
    state = node.dis0[E[1]] + arc.length.partial(cell, V.old, alpha)
  }else{
    V.old = V2
    state = node.dis0[E[2]] + arc.length.partial(cell, V.old, alpha)
  }
  cell.pt.true[i] = state
}
cell.pt.true = cbind(seq(1,n,1),cell.pt.true)
order.true = cell.pt.true[order(cell.pt.true[,2]),1]



#opt.proj1 = tree.project(X.low, opt.tree1, lambda2 = 0.01)
# initial cell
inital.cell.opt.index1 = which(apply(opt.G$Y-True.G$Y[,inital.cell.true.index], 2, sum.squares) == min(apply(opt.G$Y-True.G$Y[,inital.cell.true.index], 2, sum.squares)))

mat1 = Adjacency.mat(opt.G)
node.dis1 = spath(mat1, inital.cell.opt.index1)


# pseudotime calculation based on projection
n = dim(X)[2]
index.edge1 = proj.opt$I
cell.pt1 = rep(0,n)
for ( i in 1:n){
  e = proj.opt$E[i]
  E = index.edge1[,i]
  cell = proj.opt$fi[,i]
  V1 = opt.G$Y[,E[1]]
  V2 = opt.G$Y[,E[2]]
  alpha = opt.G$alpha[,e]
  if ( node.dis1[E[1]] < node.dis1[E[2]]){
    V.old = V1
    state = node.dis1[E[1]] + arc.length.partial(cell, V.old, alpha)
  }else{
    V.old = V2
    state = node.dis1[E[2]] + arc.length.partial(cell, V.old, alpha)
  }
  cell.pt1[i] = state
}
cell.pt1 = cbind(seq(1,n,1),cell.pt1)
order1 = cell.pt1[order(cell.pt1[,2]),1]

r1.k = cor(cell.pt1[,2], cell.pt.true[,2], method = "kendall")
r1.s = cor(cell.pt1[,2], cell.pt.true[,2], method = "spearman")
r1.p = cor(cell.pt1[,2], cell.pt.true[,2], method = "pearson")

cor.our = c(r1.k, r1.s, r1.p)
output = cor.our

### compare with linear tree
opt.proj.linear = tree.project(X, opt.tree, lambda2)
# initial cell
inital.cell.opt.index2 = which(apply(opt.tree$Y-True.G$Y[,inital.cell.true.index], 2, sum.squares) == min(apply(opt.tree$Y-True.G$Y[,inital.cell.true.index], 2, sum.squares)))

# pseudotime calculation based on projection
n = dim(X)[2]
mat2 = Adjacency.mat.linear(opt.tree$Y, opt.tree$N)
node.dis2 = spath.linear(mat2, inital.cell.opt.index2)
# pseudotime calculation based on projection
index.edge2 = opt.proj.linear$I
cell.pt.linear = rep(0,n)
for ( i in 1:n){
  E = index.edge2[,i]
  cell = opt.proj.linear$Phi[,i]
  V1 = opt.tree$Y[,E[1]]
  V2 = opt.tree$Y[,E[2]]
  if ( node.dis2[E[1]] < node.dis2[E[2]]){
    V.old = V1
    state = node.dis2[E[1]] + sqrt(sum.squares(cell-V.old))
  }else{
    V.old = V2
    state = node.dis2[E[2]] + sqrt(sum.squares(cell-V.old))
  }
  cell.pt.linear[i] = state
}
cell.pt.linear = cbind(seq(1,n,1), cell.pt.linear)
order2 = cell.pt.linear[order(cell.pt.linear[,2]),1]

r2.k = cor(cell.pt.linear[,2], cell.pt.true[,2], method = "kendall")
r2.s = cor(cell.pt.linear[,2], cell.pt.true[,2], method = "spearman")
r2.p = cor(cell.pt.linear[,2], cell.pt.true[,2], method = "pearson")

cor.our2 = c(r2.k, r2.s, r2.p)
output2 = cor.our2

