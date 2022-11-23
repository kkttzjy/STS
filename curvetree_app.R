### load dataset 
library(readxl)
raw_ct <- read_xls('mmc4.xls', 'Sheet1')


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

#::install("destiny")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

#BiocManager::install("Biobase")
library(destiny)
library(Biobase)

# as.data.frame because ExpressionSets
ct <- as.ExpressionSet(as.data.frame(raw_ct))
ct

### data cleaning
num_cells <- gsub('^(\\d+)C.*$', '\\1', ct$Cell)
ct$num_cells <- as.integer(num_cells)
# cells from 2+ cell embryos
have_duplications <- ct$num_cells > 1
# cells with values <= 28
normal_vals <- apply(exprs(ct), 2, function(smp) all(smp <= 28))
cleaned_ct <- ct[, have_duplications & normal_vals]

### normalization 
housekeepers <- c('Actb', 'Gapdh')  # houskeeper gene names

normalizations <- colMeans(exprs(cleaned_ct)[housekeepers, ])

guo_norm <- cleaned_ct
exprs(guo_norm) <- exprs(guo_norm) - normalizations

library(destiny)
#data(guo_norm)
dm <- DiffusionMap(guo_norm)
plot(dm)

palette(cube_helix(6)) # configure color palette

plot(dm,
     pch = 20,             # pch for prettier points
     col_by = 'num_cells', 
     legend_main = 'Cell stage')

X.low = t(eigenvectors(dm))

n = ncol(X.low)
p = nrow(X.low)




#### curved tree algorithm
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

#get edge list from adjancency matrix
E2N = function(E) {
  k = nrow(E)
  ind = which(E & upper.tri(E))
  N = rbind(((ind-1) %% k) + 1, floor((ind-1) / k) + 1)
  return(N)
}

## calculate the projection on the principal curve for each data point projected to the curved edge
pc = function(Y1, Y2, alpha, lambda){
  n = length(lambda)
  g = sapply(seq(1,n,1), function(x) Y1+lambda[x]*(Y2-Y1))
  h = sapply(seq(1,n,1), function(x) lambda[x]*(lambda[x]-1)*alpha)
  f = g+h
  return(f)
}


### curved tree update
tree.linear.edge.project = function(X, Y1, Y2) {
  n = dim(X)[2]
  p = dim(X)[1]
  L = pmax(0, pmin(1, colSums((X - Y1) * (Y2 - Y1)) / sum.squares(Y2 - Y1)))
  Phi = Y1 %*% t(rep(1, n)) + (Y2 - Y1) %*% t(L)
  D2 = apply(X - Phi, 2, sum.squares)
  return (list(L = L, Phi = Phi, D2 = D2))
}



tree.linear.project = function(X, G) {
  p = dim(X)[1]
  n = dim(X)[2]
  k = dim(G$Y)[2]
  I = matrix(NA, 2, n)
  L = rep(NA, n)
  Phi = matrix(NA, p, n)
  D2 = rep(Inf, n)
  for (e in 1:(k-1)) {
    i = G$N[1, e]
    j = G$N[2, e]
    proj = tree.linear.edge.project(X, G$Y[,i], G$Y[,j])
    ind = proj$D2 < D2
    I[, ind] = c(i, j)
    L[ind] = proj$L[ind]
    Phi[, ind] = proj$Phi[, ind]
    D2[ind] = proj$D2[ind]
  }
  return(list(I = I, L = L, Phi = Phi, D2 = D2))
}




tree.piece.project = function(X, a, b){
  if (is.null(dim(X)) == T){
    L = pmax(0, pmin(1, sum((X - a) * (b - a)) / sum.squares(b - a)))
    Phi = a  + (b - a) *t(L)
    D2 = sum(X - Phi)
    MAD = median(abs(X-Phi))
  } else{
    n = dim(X)[2]
    L = pmax(0, pmin(1, colSums((X - a) * (b - a)) / sum.squares(b - a)))
    Phi = a %*% t(rep(1, n)) + (b - a) %*% t(L)
    D2 = apply(X - Phi, 2, sum.squares)
    MAD = apply(abs(X-Phi), 2, median)
  }
  return(list(Phi = Phi, D2 = D2, MAD = MAD))
}

tree.curved.edge.project = function(X, Y1, Y2, alpha,lambda2) {
  n = dim(X)[2]
  p = dim(X)[1]
  #lambda = sqrt(apply(f-matrix(rep(Y1,n),p,n), 2, sum.squares)/sum.squares(Y2-Y1))
  #lambda.ord = lambda[order(lambda)]
  #knots.index =  unlist(lapply(lambda.knots, function(i) which(lambda == i)[1]))
  lambda.knots = seq(0,1,length.out = 21)
  nknots = length(lambda.knots)
  #knots.index =  unlist(lapply(lambda.knots, function(i) which(lambda == i)[1]))
  #f.knots = f[,knots.index]
  f.knots = matrix(unlist(lapply(seq(1,nknots),function(i) pc(Y1,Y2,alpha,lambda.knots[i]))), nrow=p, ncol=nknots)
  D2.all = sapply(seq(1,nknots-1,1), function(i) tree.piece.project(X, f.knots[,i], f.knots[,i+1])$D2, simplify=T)
  D2.min = apply(D2.all, 1, min)
  MAD.all = sapply(seq(1,nknots-1,1), function(i) tree.piece.project(X, f.knots[,i], f.knots[,i+1])$MAD, simplify=T)
  MAD.min = apply(MAD.all, 1, min)
  project.index = unlist(sapply(seq(1,n,1), function(i) which(D2.all[i,] == D2.min[i])[1]))
  phi.all = lapply(seq(1,nknots-1,1), function(i) tree.piece.project(X, f.knots[,i], f.knots[,i+1])$Phi)
  phi.min =  sapply(seq(1,n,1), function(i) phi.all[[project.index[i]]][,i])
  #lambda_new = sqrt(apply(phi.min-matrix(rep(Y1,n),p,n), 2, sum.squares)/sum.squares(Y2-Y1))
  lambda_new = lambda.knots[project.index]
  D2.p = lambda2*sum.squares(Y1-Y2)
  #sapply(1:n, function(i) arc.length(phi.min[,i],Y1,Y2,G$alpha))/arc.length(Y2,Y1,Y2,G$alpha)
  return (list(lambda = lambda_new, fi = phi.min, D2 = D2.min, MAD = MAD.min, D2.p = D2.p))
}

tree.curved.project = function(X, G, lambda2) {
  p = dim(X)[1]
  n = dim(X)[2]
  k = dim(G$Y)[2]
  I = matrix(NA, 2, n)
  L = rep(NA, n)
  fi = matrix(NA, p, n)
  D2 = rep(Inf, n)
  lambda = rep(Inf, n)
  SD = matrix(NA, p, k-1)
  MAD = rep(Inf, n)
  D2.p = rep(Inf, k-1)
  E = rep(NA,n)
  if (k==2){
    for (e in 1:(k-1)) {
      i = G$N[1, e]
      j = G$N[2, e]
      proj = tree.curved.edge.project(X, G$Y[,i], G$Y[,j], G$alpha, lambda2)
      ind = proj$D2 < D2
      I[, ind] = c(i, j)
      E[ind] = e
      fi[, ind] = proj$fi[, ind]
      lambda[ind] = proj$lambda[ind]
      D2[ind] = proj$D2[ind]
      MAD[ind] = proj$MAD[ind]
      D2.p[e] = proj$D2.p
    }
  } else{
    for (e in 1:(k-1)) {
      i = G$N[1, e]
      j = G$N[2, e]
      proj = tree.curved.edge.project(X, G$Y[,i], G$Y[,j], G$alpha[,e], lambda2)
      ind = proj$D2 < D2
      I[, ind] = c(i, j)
      E[ind] = e
      fi[, ind] = proj$fi[, ind]
      lambda[ind] = proj$lambda[ind]
      D2[ind] = proj$D2[ind]
      MAD[ind] = proj$MAD[ind]
      D2.p[e] = proj$D2.p
    }
  }
  return(list(I = I, lambda = lambda, fi = fi, D2 = D2, MAD = MAD, D2.p = D2.p, E = E))
}

### line search
#general numerical gradient function
grad.num = function(par, func, eps = 1e-6, ...) {
  par_new = par
  g = rep(0, length(par))
  s = func(par, ...)
  for (i in 1:length(par)) {
    par_new[i] = par[i] + eps
    s_new = func(par_new, ...)
    g[i] = (s_new - s) / eps
  }
  return (g)
}

### line search
#general numerical gradient function
grad.num = function(par, func, eps = 1e-6, ...) {
  par_new = par
  g = rep(0, length(par))
  s = func(par, ...)
  for (i in 1:length(par)) {
    par_new[i] = par[i] + eps
    s_new = func(par_new, ...)
    g[i] = (s_new - s) / eps
  }
  return (g)
}

#general line search function
line_search = function(par, func, grad = NULL, alpha_max = 1, c = 0.1, tau = 0.9, maxit = 100, diag = 0, ...) {
  s = func(par, ...)
  if (is.null(grad)) {
    g = grad.num(par, func, ...)
  } else {
    g = grad(par, ...)
  }
  if (diag > 0) {
    print("S")
    print(s)
    print("par")
    print(par)
    print("g")
    print(g)
  }
  alpha_grad = alpha_max
  par_min = par
  s_min = s
  for (it in (1:maxit)) {
    par_new = par - alpha_grad * g
    s_new = func(par_new, ...)
    if (diag > 1) print(c(s, s_new, alpha_grad))
    if (s_new < s){
      s_min = s_new
      par_min = par_new
    }
    if (s - s_new > c * alpha_grad * sum(g^2)) break
    alpha_grad = tau * alpha_grad
  }
  if (diag > 0) print(alpha_grad)
  if (it == maxit) {
    if (s_min < s){
      return(par_min)
    }
    stop("line search failed")}
  return (par_new)
}


G_score_y = function(par, X, G, lambda2, phi, i){
  G_new = G
  G_new$Y[,i] = par
  proj = tree.curved.project(X, G_new, lambda2)
  f = proj$fi
  k = dim(G$Y)[2]
  if(k == 2){
    s = sum.squares(X-f)+phi*sum(abs(G_new$alpha))+sum(proj$D2.p)
  }else{
    s = sum.squares(X-f)+phi*sum(apply(abs(G_new$alpha), 2, sum))+sum(proj$D2.p)
  }
  return(s)
}

G_score_alpha = function(par, X, G, lambda2, phi, i){
  G_new = G
  G_new$alpha[,i] = par
  proj = tree.curved.project(X, G_new, lambda2)
  f = proj$fi
  k = dim(G$Y)[2]
  if(k == 2){
    s = sum.squares(X-f)+phi*sum(abs(G_new$alpha))+sum(proj$D2.p)
  }else{
    s = sum.squares(X-f)+phi*sum(apply(abs(G_new$alpha), 2, sum))+sum(proj$D2.p)
  }
  return(s)
}

opt.score = function(X, f, phi, G){
  k = dim(G$Y)[2]
  if(k == 2){
    sum.squares(X-f)+phi*sum(abs(G$alpha))
  }else{
    sum.squares(X-f)+phi*sum(apply(abs(G$alpha), 2, sum))
  }
}

highd.plot.c = function(G, X, proj = NULL, trueplot = F, True.G = NULL){
  set.seed(12345)
  k = dim(G$Y)[2]
  n = dim(X)[2]
  PCA = prcomp(t(X), center = F, scale. = F)
  r = PCA$rotation
  #G.2d = t(t(G$Y)%*%r)
  G.2d = G$Y
  
  #plot(t(G.2d[1:2,]), col = "blue", xlim = xlimit, ylim = ylimit)
  #plot(t(G.2d[1:2,]), col = "blue")
  xlimit = range(c(G.2d[1,], X[1,]))
  ylimit = range(c(G.2d[2,], X[2,]))
  plot(X[1,],X[2,], col = "red", pch = 8,xlim = xlimit, ylim = ylimit,
       xlab =  "PC1", ylab = "PC2")
  for (e in 1:(k-1)) {
    if (k == 2){
      pcc_e = pc(G.2d[,G$N[, e][1]], G.2d[,G$N[, e][2]] , G$alpha, lambda = seq(0, 1, length.out = 100))
    }else{
      pcc_e = pc(G.2d[,G$N[, e][1]], G.2d[,G$N[, e][2]] , G$alpha[,e], lambda = seq(0, 1, length.out = 100))
    }
    xlimit = range(c(pcc_e[1,], X[1,]))
    ylimit = range(c(pcc_e[2,], X[2,]))
    lines(pcc_e[1,], pcc_e[2,], col = "blue")
    
    #plot(pc(G.2d[1, G$N[, e]], G.2d[2, G$N[, e]] , G$alpha, G$beta, lambda = seq(0, 1, length.out = 100)), type = "l", col = "black")
    #lines(t(G.2d[, G$N[, e]]), col = "black")
  }
  points(G.2d[1,], G.2d[2,], col = "blue")
  #points(PCA$x[,1],PCA$x[,2], col = "red", pch = 8)
  #plot(X[1,],X[2,], col = "red", pch = 8)
  if (trueplot == T){
    kt = dim(G.2d)[2]
    for (e in 1:(kt-1)) {
      pcc = pc(True.G$Y[,True.G$N[, e][1]], True.G$Y[,True.G$N[, e][2]] , True.G$alpha[,e], lambda = seq(0, 1, length.out = 100))
      #pcc = t(t(pcc)%*%r)
      lines(pcc[1,], pcc[2,], col = "black")
      points(True.G$Y[1,], True.G$Y[2,], col = "black")
    }
  }
  
}


## optimal curved tree
opt.tree1 = readRDS("opt_tree_line_search.rds")

G = opt.tree1
p = dim(G$Y)[1]
k = dim(G$Y)[2]

G$alpha=matrix(0.01, nrow = p, ncol = k)
X.low = t(eigenvectors(dm))[1:5,]
X = X.low
#G = True.G
highd.plot.c(G = G, X = X)

lambda2 = 0.1
c = 0.5
tau = 0.5
phi = 10
maxit = 100
tol = 1e-3
lambda1 = 1
lambda3 = 1

S_old = Inf
#highd.plot.c(G, X)
for (it in 1:maxit) {
  for (i in 1:k) {
    par2 = line_search(G$alpha[,i], G_score_alpha, NULL, c = c, tau = tau, X = X, G = G, phi = phi, lambda = lambda2, i = i,diag =2)
    G$alpha[,i] = par2
  }
  highd.plot.c(G = G, X = X)
  
  for (i in 1:k) {
    par1 = line_search(G$Y[,i], G_score_y, NULL, c = c, tau = tau, X = X, G = G, phi = phi, lambda = lambda2, i = i,diag =2)
    G$Y[,i] = par1
  }
  highd.plot.c(G = G, X = X)
  
  proj = tree.curved.project(X, G, lambda2)
  S = opt.score(X, proj$fi, phi, G)+sum(proj$D2.p)
  if (S > S_old) stop("S > S_old")
  if (S_old - S < tol) break
  S_old = S
  proj = tree.curved.project(X, G, lambda2)
}


pdf("linesearch.curved.pdf")
highd.plot.c(G = G, X = X)
dev.off()

save(G, file="linesearch_opt.G.Rdata")


load("linesearch_opt.G.Rdata")
opt.G = G
proj.opt = tree.curved.project(X, opt.G, lambda2)
highd.plot.c(G = opt.G, X = X)

library(ggplot2)
pt = cell.stage.num
plotdata = data.frame(cbind(t(X.low), pt))

k = dim(opt.G$Y)[2]

for (e in 1:(k-1)) {
  if (k == 2){
    pcc_e = pc(opt.G$Y[,opt.G$N[, e][1]], opt.G$Y[,opt.G$N[, e][2]] , opt.G$alpha, lambda = seq(0, 1, length.out = 100))
  }else{
    pcc_e = pc(opt.G$Y[,opt.G$N[, e][1]], opt.G$Y[,opt.G$N[, e][2]] , opt.G$alpha[,e], lambda = seq(0, 1, length.out = 100))
  }
  pcc_e = data.frame(t(pcc_e))
  if (e == 1){
    sp = ggplot(data = pcc_e, aes(x = DC1, y = DC2))
  }
  # TY = data.frame(x1 = opt.G$Y[1, opt.G$N[, e][1]], 
  #                 x2 = opt.G$Y[1, opt.G$N[, e][2]], 
  #                 y1 = opt.G$Y[2, opt.G$N[, e][1]], 
  #                 y2 = opt.G$Y[2, opt.G$N[, e][2]])
  TP = data.frame(t(opt.G$Y[1:2, opt.G$N[, e]]))
  sp = sp + geom_point(aes(x = DC1, y = DC2), color = "#333333", data = TP) +
    geom_line(data = pcc_e, aes(x = DC1, y = DC2), size = 1)
}

sp = sp + geom_point(data = plotdata, aes(x = DC1, y = DC2,  color = pt))
sp = sp + scale_color_gradient(low="blue", high="red")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(sp)

pdf("curved.linesearch.plot.pdf")
highd.plot.c(opt.G, X.low)
print(sp)
dev.off()

# initial cell
### cell stage
cell.stage = dm$Cell
cell.stage = unlist(strsplit(cell.stage, "C"))[seq(1,2*length(cell.stage),2)]
cell.stage.num = as.numeric(cell.stage)



start.cell = which(cell.stage == "2")
start.edge = proj.opt$I[,start.cell][,1]
start.edge.dis = sqrt(apply(opt.G$Y[,start.edge]-apply(X.low[,start.cell],1,mean), 2, sum.squares))
initial.cell.opt.index1 = start.edge[which(start.edge.dis == min(start.edge.dis))]
initial.cell.opt.pt1 = min(start.edge.dis)
Tree.temp = opt.G
k = ncol(opt.G$Y)
Tree.temp$Y = cbind(opt.G$Y, proj.opt$Phi[,start.cell])
mat1 = Adjacency.mat(opt.G)
node.dis1 = spath(mat1, initial.cell.opt.index1)

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
  cell.pt1[i] = state + initial.cell.opt.pt1
}
cell.pt1[start.cell] = initial.cell.opt.pt1
cell.pt1 = cbind(seq(1,n,1), cell.pt1)
order1 = cell.pt1[order(cell.pt1[,2]),1]


r1.k = cor(cell.pt1[,2], cell.stage.num, method = "kendall")
r1.s = cor(cell.pt1[,2], cell.stage.num, method = "spearman")
r1.p = cor(cell.pt1[,2], cell.stage.num, method = "pearson")

cor.our = c(r1.k, r1.s, r1.p)
output = cor.our

cell.pt = cell.pt1[,2]
cell.pt = as.character(cell.pt)
cell.pt = as.numeric(cell.pt)
plotdata = data.frame(cbind(cell.stage.num, cell.pt))
colnames(plotdata) = c("true", "estimated")
library(ggplot2)
pdf("curved.plot.our.pt.pdf")
p1 = ggplot(plotdata, aes(x = true, y = estimated)) +
  geom_point()+
  labs(
    x = "True time",
    y = "Estimated pseudotime"
  )
print(p1)
dev.off()

