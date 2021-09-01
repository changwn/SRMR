
simu_outlier <- function(vv, beta=c(1.5, -1.5), distOut = 0.5){
  out = rep(0, length(vv))
  i = 1
  while(i<= length(vv)){
    tmp = sample(seq(-8, 8,0.01),1)
    if(abs(tmp - beta[1]*vv[i])>distOut & abs(tmp - beta[2]*vv[i])>distOut ){ #reject if tmp close to mixture lines.
      out[i] = tmp
      i = i+1
    }
  }
  return(out)
}

library(MASS)
simulateData_outlier <- function(n=c(40,40,20), beta=c(1.5,-1.5), coordinate=c(1,1,-1,-1)){
  cl = c(rep(1,35), rep(4, 5), 
         rep(2,35), rep(5, 5), 
         rep(3, 20))
  
  mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
  colnames(mat) = c('y','x')
  loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:100
  mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,0.1)
  mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,0.1)
  # mat[loca3,1] = abs(beta[1] * mat[loca3,2]) + abs(beta[2] * mat[loca3,2]) + rnorm(20,0,3)
  mat[loca3,1] = simu_outlier(mat[loca3, 2], beta, distOut=2)
  mat = mat[,2:1] #change order
  
  # spatial coordinate
  ccc = mvrnorm(n[1], mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
  ccc2 = mvrnorm(n[2], mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
  ccc.o1 = mvrnorm(n[3]/2, mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
  ccc.o2 = mvrnorm(n[3]/2, mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
  ccc.n1 = rbind(ccc[1:35,], ccc2[36:40,])
  ccc.n2 = rbind(ccc2[1:35,], ccc[36:40,])
  xy = rbind(ccc.n1, ccc.n2,ccc.o1, ccc.o2)
  
  return(list(mat=mat, xy=xy, cl=cl))
}



#debug default parameter
# formula= as.formula("y~x")
# # data
# nit=20
# nc=2
# rlr_method="ltsReg"
# Cdn=xy
# lamb=5
# showPlot=F
RobSpaReg<- function(formula,data, nit=20,nc=2,rlr_method="ltsReg", Cdn=xy, lamb=5, showPlot=F){
  mycall = match.call();
  res_list=vector("list",nit)
  ooo_list=vector("list",nit)
  nx=ncol(data)-1; nobs = nrow(data)            # number of observations
  for(jj in 1:nit){
    flag=0; outliers=c(); ccc=0
    www=NULL
    for(ii in 1:nc){
      inds_random=sample(1:nobs,((nx+1)*2*5))
      if(rlr_method =="lmRob"){
        mod_tmp = try(lmRob(formula=formula, data=data[inds_random,],control = lmRob.control(weight=c("Bisquare", "Bisquare"))),silent=TRUE)
      }
      if(rlr_method =="lmrob"){
        mod_tmp = try(lmrob(formula=formula, data=data[inds_random,]),silent=TRUE)
      }
      if(rlr_method =="ltsReg"){
        mod_tmp = try(ltsReg(formula=formula, data=data[inds_random,]),silent=TRUE)
      }
      #www=cbind(www,dnorm(predict(mod_tmp, newdata=data),mean=0,sd=summary(mod_tmp)$sigma))
      www=cbind(www,dnorm(abs(cbind(1,as.matrix(data)) %*% matrix(c(mod_tmp$coefficients,-1),ncol=1))[,1],mean=0,sd=summary(mod_tmp)$sigma))
    }
    nc_tmp=nc
    while(flag==0 ){
      inds_in=1:nobs; outliers.old=outliers
      ccc=ccc+1
      outlier_list1=outlier_list=vector("list", nc_tmp)
      for(j in 1:nc_tmp){
        inds_in_tmp=which(apply(www, 1, which.max)==j)
        #ltsres = lm(formula=formula, data=data[inds_in_tmp,])
        if(rlr_method =="lmRob"){
          ltsres = try(lmRob(formula=formula, data=data[inds_in_tmp,],control = lmRob.control(weight=c("Bisquare", "Bisquare"))),silent=TRUE)
        }
        if(rlr_method =="lmrob"){
          ltsres = try(lmrob(formula=formula, data=data[inds_in_tmp,]),silent=TRUE)
        }
        if(rlr_method =="ltsReg"){
          ltsres = try(ltsReg(formula=formula, data=data[inds_in_tmp,]),silent=TRUE)
        }
        if(!inherits(ltsres, "try-error")){
          res12 = abs(cbind(1,as.matrix(data)) %*% matrix(c(ltsres$coefficients,-1),ncol=1))[,1] ##could we change it to nx!!!.
          oo1=inds_in_tmp[DeOut(res12[inds_in_tmp],"hampel")]
          if(rlr_method=="lmRob"){
            #oo2=inds_in_tmp[which(ltsres$T.M.weights==0)]
            oo2=inds_in_tmp[which(ltsres$M.weights==0)]
          }
          if(rlr_method=="ltsReg"){
            oo2=inds_in_tmp[which(ltsres$lts.wt==0)]
          }
          if(rlr_method=="lmrob"){
            oo2=inds_in_tmp[which(ltsres$rweights==0)]
          }
          #outlier_list[[j]] = oo1##i
          #outlier_list[[j]] = oo2##ii
          outlier_list[[j]] = intersect(oo1,oo2)##iii
          #outlier_list[[j]] = unique(c(oo1,oo2))##c
        }
      } #end-for_nc
      # print(outlier_list)
      outliers=Reduce(c,outlier_list)
      if(length(outliers)>0){inds_in=c(1:nobs)[-outliers]}
      if (showPlot==T){
        ccol = rep(1, nobs); ccol[outliers] = 2; plot(data, col=ccol,pch=16) #debug
      }
      # fres = flexmix_2(formula,data1=data[inds_in,],k=nc,mprior=0.1)            # mix_reg. Output: posterior, logLik, 
      # fres@cluster
      SpaRes = SpatialRegKmeans(dat=data[inds_in,], ncl=nc, iter.max = 100L, epsilon=1e-4, Cdn_f=Cdn[inds_in,], verbose=F, lambda=lamb, showPlot=F, rob=T, inside=inds_in)
      # SpaRes.only = SpatialRegKmeans(dat=data[inds_in,], ncl=nc, iter.max = 100L, epsilon=1e-4, Cdn=xy[inds_in,], verbose=F, lambda=0, showPlot=F)
      # SpaRes$clusterMem
      # print(rbind(inds_in, SpaRes$clusterMem))
      # ### backup step for catching spatial outliers.
      # print(rbind(inds_in, SpaRes$clusterMem, SpaRes.only$clusterMem))
      # spa.outlier = xxx(rbind(SpaRes$clusterMem, SpaRes.only$clusterMem))
      if (showPlot==T){
        # print(SpaRes$outlier_spa)
        # ccol = rep(1, nobs); ccol[outliers] = 2; ccol[SpaRes$outlier_spa]=3; plot(data, col=ccol,pch=16) #debug
      }
      # print(outliers)
      # print(SpaRes$outlier_spa)
      outliers = unique(sort(c(outliers, SpaRes$outlier_spa)))
      # print(outliers)
      # print('-----------')

      # nc_tmp=fres@k
      # www = posterior(fres, newdata=data, unscale=TRUE)
      www = calPosterior(SpaRes, newdata=data, newCdn=Cdn, ncl=nc, lambda=lamb)
      # pwww = t(www); colnames(pwww) = c(1:ncol(pwww)) ; barplot(pwww, beside=TRUE, col=c("lightblue","red"),cex.names=0.5, main='www: posterior')
      if(ccc>10  ){flag=1}
      if(length(outliers)==length(outliers.old)& ccc>1){
        if(length(outliers)==sum(outliers==outliers.old)){
          flag=1;
        }
      }
    }#end-while_flag
    #print(sort(unique(outliers)))
    #if(fres@k == nc & ccc<11){
    if(nrow(SpaRes$centroid) == nc & ccc<11){
      # res_list[[jj]]=fres # mix_reg result. 
      res_list[[jj]]=SpaRes
      ooo_list[[jj]]=sort(unique(outliers))
      #print(dim(SpaRes$hy_posterior))
      # print(length(sort(unique(outliers))))
      #print(length(sort(outliers)))
    }
    #print(paste("The number of iteraction is", ccc))
  } #end-jj
  ooo_list=ooo_list[][which(sapply(res_list, length)>0)]
  res_list=res_list[][which(sapply(res_list, length)>0)]
  #llik=sapply(res_list, function(x)x@logLik)
  llik=sapply(res_list, function(x)x$loss)
  res_list=res_list[][order(llik,decreasing=TRUE)]
  ooo_list=ooo_list[][order(llik,decreasing=TRUE)]
  
  list1=sapply(ooo_list, function(x)list(1*((1:nobs)%in%x)))
  opt.ind=which.min(sapply(list1, function(x)sum((x-Reduce("+",list1)/length(list1))^2)))
  opt.fres=res_list[[opt.ind]]
  
  outliers=ooo_list[[opt.ind]]
  inds_in=setdiff(1:nobs,outliers)
  #coffs_final=rbind(parameters(opt.fres),apply(opt.fres@posterior$scaled,2,sum)/length(inds_in))
  coffs_final=rbind(opt.fres$alpha, opt.fres$beta); rownames(coffs_final) = c('coef.(Intercept)','coef.x'); 
  coffs_final=rbind(coffs_final, apply(opt.fres$hy_posterior,2,sum)/length(inds_in))
  cates=vector("numeric",nobs)
  cates[inds_in]=opt.fres$clusterMem #opt.fres@cluster
  cates[outliers]=-1
  
  wij=matrix(NA,nrow=nobs,ncol=nc)
  wij[inds_in,]=opt.fres$hy_posterior #opt.fres@posterior$unscaled
  
  #pvals_final=sapply(refit(opt.fres)@components[[1]], function(x)(x[-1,4]))
  # result = new("RobMixReg", inds_in=inds_in,indout=outliers,ctleclusters=cates,compcoef=coffs_final,comppvals=pvals_final,compwww=wij)
  result = new("RobMixReg", inds_in=inds_in,indout=outliers,ctleclusters=cates,compcoef=coffs_final,compwww=wij)
  result@call <- mycall
  result
}

calPosterior <- function(SpaRes, newdata=data, newCdn=xy, ncl=nc, lambda=1){
  mmm = matrix(0, nrow(newdata), ncl)
  #y = newdata$y; x = newdata$x
  y = newdata[,1]; x = newdata[,2]  # assume 1st col is y, 2nd col is x.
  alpha = SpaRes$alpha
  beta = SpaRes$beta
  s.center = SpaRes$centroid
  
  ##### update C_k
  Y = matrix(y, length(y), ncl)
  A = matrix(1, length(y), 1) %*% alpha
  X = matrix(x, length(x), 1) %*% beta
  rownames(Y) = rownames(X) = names(y)
  #Z = (Y-X)^2
  # clusterMem = apply((Y - A - X)^2, 1, which.min) #residual
  resi = abs(Y - A - X)
  resi_scale = 1 - resi / rowSums(resi) #sum is zero?
  
  D = matrix(0, length(y), ncl)
  for(k in 1:ncl){
    center.tmp = s.center[k,]
    D[, k] = apply(newCdn, 1, function(x){norm(matrix(x - center.tmp, nrow=1), 'f')})
  }
  D_scale = 1 - D / rowSums(D)
  
  # hybrid: Ci = p(Ci | xi, beta) + lambda * p(Ci | xi, center)
  p_hybrid = resi_scale + lambda * D_scale
  
  return(p_hybrid)
}

# spatial regression-wised kmeans.
# x : n by 2 matrix, only support 2-dim features.
# ncl: number of clusters.
# iter.max: maximum iteration number, default is 100.
# epsilon: converged threshold, default is 1e-4.
# Cdn: coordinate of data, which is a n by 2 matrix, n is number of observations.
# verbose: if ture, will print iteration result.
# lambda: hyperparameter of hybridzation, which is to balance regression-based probability and spatial-based probability.
SpatialRegKmeans <- function(dat, ncl, iter.max = 100L, epsilon=1e-4, Cdn_f=NULL, verbose=T, lambda=1, showPlot=F, rob=F, inside)
{
  
  n_sample = nrow(dat) ; 
  y = dat[,2]; x = dat[,1] #order is important: 1st column is x, 2nd column is y.
  
  if(ncol(Cdn_f) != 2) stop('The coordinate should be 2-dim.')
  #clusterMem = rep(1, n_sample) # give random cluster
  #clusterMem = sapply(clusterMem, function(x){sample(1:ncl, 1, replace = T)})
  # d_mat <- as.dist(W)
  # hclust.res <- hclust(d_mat, method = "complete" )
  # clusterMem <- cutree(hclust.res, k=ncl) # cut tree into 5 clusters
  kmeans.res <- kmeans(Cdn_f, ncl)
  clusterMem <- kmeans.res$cluster
  km.center <- kmeans.res$centers
  if(verbose==T) {print('Init random cluster membership.'); print(table(clusterMem))}
  
  if(showPlot==T) plot(dat, col=clusterMem, main='Random Init') 
  # hist(W, main='Dist of distance matrix')
  if(showPlot==T) plot(Cdn_f, main='Coordinate')
  
  
  loss = Inf
  loss_small = Inf
  clusterMem_good = clusterMem
  beta_good = matrix(0, 1, ncl)
  count = 0
  beta = matrix(0, 1, ncl)
  alpha = matrix(0, 1, ncl)
  spa_outlier = list()
  # while(loss > epsilon | count < iter.max){
  while(count < iter.max){
    loss = 0
    count = count + 1
    if(verbose == T) print(paste('iter=', count, sep='') )
    
    if(length(table(clusterMem)) < ncl | min(table(clusterMem)) == 1){
      if(verbose == T) print('Reshuffled.')
      clusterMem = sapply(clusterMem, function(x){sample(1:ncl, 1, replace = T)}) # reshuffle 
    }
    
    for(k in 1:ncl){ #update spatial center
      km.center[k, ] = apply(Cdn_f[clusterMem==k,], 2, mean)
    }
    
    # update B_k
    for(k in 1:ncl){
      if(verbose==T)print(paste('current k:',k, sep=' '))
      candi = which(clusterMem == k)
      
      #method 4: optim package
      n_slot = choose(length(candi),2)
      # W_slot = matrix(0, n_slot, 1)
      y_slot = matrix(0, n_slot, 1)
      x_slot = matrix(0, n_slot, 1)
      L = 0
      ks = 0
      for(i in 1:(length(candi)-1)){ 
        for(j in (i+1):length(candi)){
          id_i = candi[i]
          id_j = candi[j]
          # L = L + W[id_i, id_j] * (y[id_i]-y[id_j] + beta_k*(x[id_j]-x[id_i]) )^2
          ks = ks + 1
          # W_slot[ks] = W[id_i, id_j] #^ pow # increase importance of distance.
          y_slot[ks] = y[id_i] - y[id_j]
          x_slot[ks] = x[id_j] - x[id_i]
        }
      }
      lossFunc <- function(x_slot, y_slot, par){
        # J <- apply(cbind(y_slot, x_slot, W_slot), 1, function(x){(x[1] + par * x[2])^2 * x[3]} )
        # J <- apply(cbind(y_slot, x_slot, W_slot), 1, function(x){(x[1] + par * x[2])^2 + 100*x[3]} )
        J <- apply(cbind(y_slot, x_slot), 1, function(x){(x[1] + par * x[2])^2 } )
        # print(sum(apply(cbind(y_slot, x_slot, W_slot), 1, function(x){(x[1] + par * x[2])^2 } )))
        # print('----')
        # print(sum(apply(cbind(y_slot, x_slot, W_slot), 1, function(x){x[3]} ) ))
        J <- sum(J)
        return(J)
      }
      res <- optim(par = rep(0,1), fn = lossFunc, y_slot=y_slot, x_slot=x_slot, method = "BFGS")
      # res$par
      # res$value
      beta[k] = res$par
      loss = loss + res$value
    } #end-for
    if(verbose==T){print(paste('iter ', count, ', beta:', sep=''));print(beta)}
    if(verbose==T){print(paste('iter ', count, ', loss:', sep=''));print(loss)}
    
    # estimate intercept a_k
    for(k in 1:ncl){
      candi = which(clusterMem == k)
      # median(y[candi] - x[candi] * beta[k])
      # plot(y[candi] - x[candi] * beta[k])
      # boxplot(y[candi] - x[candi] * beta[k])
      alpha[k] = median(y[candi] - x[candi] * beta[k])
    }
    if(verbose==T){print(paste('iter ', count, ', alpha:', sep=''));print(alpha)}
    
    ##### update C_k
    Y = matrix(y, length(y), ncl)
    A = matrix(1, length(y), 1) %*% alpha
    X = matrix(x, length(x), 1) %*% beta
    rownames(Y) = rownames(X) = names(y)
    #Z = (Y-X)^2
    # clusterMem = apply((Y - A - X)^2, 1, which.min) #residual
    resi = abs(Y - A - X)
    resi_scale = 1 - resi / rowSums(resi) #sum is zero?
    clust_resi = apply(resi_scale, 1, which.max)
    
    D = matrix(0, length(y), ncl)
    for(k in 1:ncl){
      center.tmp = km.center[k,]
      D[, k] = apply(Cdn_f, 1, function(x){norm(matrix(x - center.tmp, nrow=1), 'f')})
    }
    D_scale = 1 - D / rowSums(D)
    cluster_spa = apply(D_scale, 1, which.max)
    # print(rbind(inds_in, clust_resi, cluster_spa))
    spa_outlier[[length(spa_outlier)+1]] = which(clust_resi != cluster_spa)
    
    # hybrid: Ci = p(Ci | xi, beta) + lambda * p(Ci | xi, center)
    p_hybrid = resi_scale + lambda * D_scale
    clusterMem = apply(p_hybrid, 1, which.max) 
    
    if(verbose==T) print(table(clusterMem))
    
    if(showPlot==T) plot(dat, col=clusterMem, main=paste('loss=',loss,sep=''))
    
    if(loss < loss_small){
      clusterMem_good = clusterMem
      loss_small = loss
      beta_good = beta
    }else{
      break;
      # cause diverse...
      # reshuffle partial based on best solution
      # n_partial = 0.01*n_sample
      # loca = sample(n_sample, n_partial)
      # clusterMem = clusterMem_good
      # clusterMem[loca] = sapply(clusterMem[loca], function(x){sample(1:ncl, 1, replace = T)})
    }
    
  } #end-while
  clusterMem_rt = clusterMem_good
  spa_outlier_c = Reduce(intersect, spa_outlier)
  spa_outlier_trans = inside[spa_outlier_c]
  
  if(length(spa_outlier_c) == 0){ #perfect situation
    return(list(clusterMem=clusterMem_rt, loss=loss_small, beta=beta, alpha=alpha, centroid=km.center, 
                hy_posterior=p_hybrid, outlier_spa=NULL))
  }
  
  # return(list(clusterMem=clusterMem_rt, loss=loss_small, beta=beta, alpha=alpha, centroid=km.center, hy_posterior=p_hybrid, outlier_spa=spa_outlier_c))
  return(list(clusterMem=clusterMem_rt[-spa_outlier_c], loss=loss_small, beta=beta, alpha=alpha, centroid=km.center, 
              hy_posterior=p_hybrid[-spa_outlier_c,], outlier_spa=spa_outlier_trans))
}