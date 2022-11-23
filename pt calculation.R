lambda2 = 0.01
c = 0.1
tau = 0.5
phi = 0.01
maxit = 200
tol = 1e-3
lambda1 = 1
lambda3 = 1
results.k = NULL
results.s = NULL
results.p = NULL
results.ss= NULL
results.mse = NULL
#id.index = c(1,6,8,9,12,13,14,15,16,18,19,20,21,23,24,25,26,27,28,29,30,32,33,35,36,37,
#             38,39,41,43,44,45,47,48,49,51,52,53,54,56,59,62,63,65,68,69,70,71,73,75,
#             76,77,78,79,80,81,84,85,86,91,93,94,95,96,97,98)
for ( id in 1:200){
  tryCatch({
  ## optimal curved tree
  load(paste0(id, "_data.Rdata"))
  load(paste0(id, "_trueG.Rdata"))
  load(paste0(id, "_opttree.Rdata"))
  # G = opt.tree1
  # p = dim(G$Y)[1]
  # k = dim(G$Y)[2]
  # 
  # G$alpha=matrix(0.5, nrow = p, ncol = k)
  # #G = True.G
  # highd.plot.c(G = G, X = X, trueplot = T, True.G = True.G)
  # 
  # 
  # S_old = Inf
  # #highd.plot.c(G, X)
  # for (it in 1:maxit) {
  #   for (i in 1:k) {
  #     par2 = line_search(G$alpha[,i], G_score_alpha, NULL, c = c, tau = tau, X = X, G = G, phi = phi, lambda = lambda2, i = i,diag =2)
  #     G$alpha[,i] = par2
  #   }
  #   highd.plot.c(G = G, X = X, trueplot = T, True.G = True.G)
  #   
  #   for (i in 1:k) {
  #     par1 = line_search(G$Y[,i], G_score_y, NULL, c = c, tau = tau, X = X, G = G, phi = phi, lambda = lambda2, i = i,diag =2)
  #     G$Y[,i] = par1
  #   }
  #   highd.plot.c(G = G, X = X, trueplot = T, True.G = True.G)
  #   
  #   proj = tree.curved.project(X, G, lambda2)
  #   S = opt.score(X, proj$fi, phi, G)+sum(proj$D2.p)
  #   if (S > S_old) stop("S > S_old")
  #   if (S_old - S < tol) break
  #   S_old = S
  #   proj = tree.curved.project(X, G, lambda2)
  # }
  # pdf(paste0(id, ".curved.pdf")) 
  # highd.plot.c(G = G, X = X, trueplot = T, True.G = True.G)
  # dev.off()
  # save(G, file=paste0(id, "_opt.G.Rdata"))
  # 
  # 
  load(paste0(id, "_opt.G.Rdata"))
  ### calculate psuedotime
  opt.G = G
  proj.opt = tree.curved.project(X, opt.G, lambda2)
  k = dim(opt.G$Y)[2]
  p = dim(opt.G$Y)[1]
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
  opt.proj.linear = tree.project(X, opt.tree1, lambda2)
  # initial cell
  inital.cell.opt.index2 = which(apply(opt.tree1$Y-True.G$Y[,inital.cell.true.index], 2, sum.squares) == min(apply(opt.tree1$Y-True.G$Y[,inital.cell.true.index], 2, sum.squares)))
  
  # pseudotime calculation based on projection
  n = dim(X)[2]
  mat2 = Adjacency.mat.linear(opt.tree1$Y, opt.tree1$N)
  node.dis2 = spath.linear(mat2, inital.cell.opt.index2)
  # pseudotime calculation based on projection
  index.edge2 = opt.proj.linear$I
  cell.pt.linear = rep(0,n)
  for ( i in 1:n){
    E = index.edge2[,i]
    cell = opt.proj.linear$Phi[,i]
    V1 = opt.tree1$Y[,E[1]]
    V2 = opt.tree1$Y[,E[2]]
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
  results = cbind(output,output2)
  r.k = c(r1.k,r2.k)
  r.s = c(r1.s,r2.s)
  r.p = c(r1.p,r2.p)
  
  ss.residual.curve = sum(apply(true.proj$fi - proj.opt$fi, 2, sum.squares))
  ss.residual.linear = sum(apply(true.proj$fi - opt.proj.linear$Phi, 2, sum.squares))
  mse.linear = sqrt(ss.residual.linear/100)
  mse.curve = sqrt(ss.residual.curve/100)

  ss.all = c(ss.residual.curve, ss.residual.linear)
  mse.all = c(mse.curve, mse.linear)
  
  results.k = rbind(results.k, r.k)
  results.s = rbind(results.s, r.s)
  results.p = rbind(results.p, r.p)
  results.ss = rbind(results.ss, ss.all)
  results.mse = rbind(results.mse, mse.all)
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
Nsim = dim(results.k)[1]

apply(results.k,2,mean)
apply(results.k,2,sd)

apply(results.p,2,mean)
apply(results.p,2,sd)

apply(results.s,2,mean)
apply(results.s,2,sd)

apply(results.ss,2,mean)
apply(results.ss,2,sd)

apply(results.mse,2,mean)
apply(results.mse,2,sd)

