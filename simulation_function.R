simulateData_outlier_sig <- function(n=c(40,40,20), beta=c(1.5,-1.5), sigma=0.05, 
                                     coordinate=c(1,1,-1,-1) 
)
{
  cl = c(rep(1,35), rep(4, 5), 
         rep(2,35), rep(5, 5), 
         rep(3, 20))
  cl_new = c(rep(1,35), rep(-1, 5), 
             rep(2,35), rep(-2, 5), 
             rep(-3, 20))
  
  mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
  colnames(mat) = c('y','x')
  loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:100
  mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
  mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
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
  
  return(list(mat=mat, xy=xy, cl=cl, cl_new=cl_new, beta=beta))
}


simulateData_outlier_n <- function(n=c(40,40,20), beta=c(1.5,-1.5), sigma=0.05, 
                                   coordinate=c(1,1,-1,-1) 
                                  )
{
  N = sum(n)
  cl = c(rep(1,N*0.35), rep(4, N*0.05), 
         rep(2,N*0.35), rep(5, N*0.05), 
         rep(3, N* 0.2))
  cl_new = c(rep(1,N*0.35), rep(-1, N*0.05), 
             rep(2,N*0.35), rep(-2, N*0.05), 
             rep(-3, N* 0.2))
  
  mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
  colnames(mat) = c('y','x')
  loca1 = which(cl == 1 | cl == 4)  
  loca2 = which(cl == 2 | cl == 5)  
  loca3 = which(cl == 3) 
  mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
  mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
  mat[loca3,1] = simu_outlier(mat[loca3, 2], beta, distOut=2)
  mat = mat[,2:1] #change order
  
  # spatial coordinate
  ccc = mvrnorm(n[1], mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
  ccc2 = mvrnorm(n[2], mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
  ccc.o1 = mvrnorm(n[3]/2, mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
  ccc.o2 = mvrnorm(n[3]/2, mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
  loca.in = which(cl == 1); loca.out = which(cl == 4)
  ccc.n1 = rbind(ccc[loca.in,], ccc2[loca.out,])
  ccc.n2 = rbind(ccc2[loca.in,], ccc[loca.out,])
  xy = rbind(ccc.n1, ccc.n2,ccc.o1, ccc.o2)
  
  return(list(mat=mat, xy=xy, cl=cl, cl_new=cl_new, beta=beta))
}


simu_outlier_k <- function(vv, beta=NULL, distOut = 0.5){
  out = rep(0, length(vv))
  i = 1
  n_comp = length(beta)
  
  while(i<= length(vv)){
    acept=TRUE
    tmp = sample(seq(-8, 8,0.01),1)
    for(j in 1:n_comp){
      if(abs(tmp - beta[j]*vv[i])<=distOut) acept=FALSE
    }
    if(acept){ #reject if tmp close to mixture lines.
      out[i] = tmp
      i = i+1
    }
  }
  return(out)
}

simulateData_outlier_k <- function(k = 2)
{
  if(! k%in% c(2,3,4) ) stop('k can only be 2,3,4!')
  if(k == 2){
    n=c(40,40,20);
    beta=c(1.5,1.0)
    sigma=0.05 
    coordinate=c(1,1,-1,-1)
    
    cl = c(rep(1,35), rep(4, 5), 
           rep(2,35), rep(5, 5), 
           rep(3, 20))
    cl_new = c(rep(1,35), rep(-1, 5), 
               rep(2,35), rep(-2, 5), 
               rep(-3, 20))
    
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:100
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
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
  }
  
  if(k == 3){
    n=c(40,40,40,20);
    beta=c(3,1.5,-3)
    sigma=0.05 
    coordinate=c(1,1, -1,-1, -1,1)
    
    cl = c(rep(1,35), rep(5, 5), 
           rep(2,35), rep(6, 5), 
           rep(3,35), rep(7, 5), 
           rep(4, 20))
    cl_new = c(rep(1,35), rep(-5, 5), 
           rep(2,35), rep(-6, 5), 
           rep(3,35), rep(-7, 5), 
           rep(-4, 20))
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:120; loca4 = 121:140
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
    mat[loca3,1] = beta[3] * mat[loca3,2] + rnorm(40,0,sigma)
    mat[loca4,1] = simu_outlier_k(mat[loca4, 2], beta, distOut=2)
    mat = mat[,2:1] #change order
    
    # spatial coordinate
    ccc = mvrnorm(n[1], mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc2 = mvrnorm(n[2], mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
    ccc3 = mvrnorm(n[3], mu = c(coordinate[5],coordinate[6]), Sigma=diag(0.1, 2,2))
    ccc.o1 = mvrnorm(7, mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc.o2 = mvrnorm(7, mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
    ccc.o3 = mvrnorm(6, mu = c(coordinate[5],coordinate[6]), Sigma=diag(0.1, 2,2))
    ccc.n1 = rbind(ccc[1:35,], ccc2[36:40,])
    ccc.n2 = rbind(ccc2[1:35,], ccc3[36:40,])
    ccc.n3 = rbind(ccc3[1:35,], ccc[36:40,])
    xy = rbind(ccc.n1,ccc.n2,ccc.n3,  ccc.o1,ccc.o2,ccc.o3)
  }
  
  if(k == 4){
    n=c(40,40,40,40,20);
    beta=c(3,1.5,-1.5,-3)
    sigma=0.05 
    coordinate=c(1,1, -1,-1, -1,1, 1,-1)
    
    cl = c(rep(1,35), rep(6, 5), 
           rep(2,35), rep(7, 5), 
           rep(3,35), rep(8, 5), 
           rep(4,35), rep(9, 5), 
           rep(5, 20))
    cl_new = c(rep(1,35), rep(-6, 5), 
           rep(2,35), rep(-7, 5), 
           rep(3,35), rep(-8, 5), 
           rep(4,35), rep(-9, 5), 
           rep(-5, 20))
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:120; loca4=121:160; loca5 = 161:180
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
    mat[loca3,1] = beta[3] * mat[loca3,2] + rnorm(40,0,sigma)
    mat[loca4,1] = beta[4] * mat[loca4,2] + rnorm(40,0,sigma)
    mat[loca5,1] = simu_outlier_k(mat[loca5, 2], beta, distOut=2)
    mat = mat[,2:1] #change order
    
    # spatial coordinate
    ccc = mvrnorm(n[1], mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc2 = mvrnorm(n[2], mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
    ccc3 = mvrnorm(n[3], mu = c(coordinate[5],coordinate[6]), Sigma=diag(0.1, 2,2))
    ccc4 = mvrnorm(n[4], mu = c(coordinate[7],coordinate[8]), Sigma=diag(0.1, 2,2))
    ccc.o1 = mvrnorm(5, mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc.o2 = mvrnorm(5, mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
    ccc.o3 = mvrnorm(5, mu = c(coordinate[5],coordinate[6]), Sigma=diag(0.1, 2,2))
    ccc.o4 = mvrnorm(5, mu = c(coordinate[7],coordinate[8]), Sigma=diag(0.1, 2,2))
    ccc.n1 = rbind(ccc[1:35,], ccc2[36:40,])
    ccc.n2 = rbind(ccc2[1:35,], ccc3[36:40,])
    ccc.n3 = rbind(ccc3[1:35,], ccc4[36:40,])
    ccc.n4 = rbind(ccc4[1:35,], ccc[36:40,])
    xy = rbind(ccc.n1,ccc.n2,ccc.n3,ccc.n4,  ccc.o1,ccc.o2,ccc.o3,ccc.o4)
  }
  
  
  return(list(mat=mat, xy=xy, cl=cl, cl_new=cl_new, beta=beta))
}


simulateData_outlier_bal <- function(n=c(40,40,20), beta=c(1.5,-1.5), sigma=0.05, 
                                     coordinate=c(1,1,-1,-1) 
)
{
  if(n[1] == 40){
    cl = c(rep(1,35), rep(4, 5), 
           rep(2,35), rep(5, 5), 
           rep(3, 20))
    cl_new = c(rep(1,35), rep(-4, 5), 
           rep(2,35), rep(-5, 5), 
           rep(-3, 20))
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:100
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
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
  }
  
  if(n[1] == 30){ 
    cl = c(rep(1,25), rep(4, 5), 
           rep(2,45), rep(5, 5), 
           rep(3, 20))
    cl_new = c(rep(1,25), rep(-4, 5), 
           rep(2,45), rep(-5, 5), 
           rep(-3, 20))
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:30;  loca2 = 31:80;  loca3 = 81:100
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(30,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(50,0,sigma)
    mat[loca3,1] = simu_outlier(mat[loca3, 2], beta, distOut=2)
    mat = mat[,2:1] #change order
    
    # spatial coordinate
    ccc = mvrnorm(n[1], mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc2 = mvrnorm(n[2], mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
    ccc.o1 = mvrnorm(5, mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc.o2 = mvrnorm(15, mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2)) #15+5
    ccc.n1 = rbind(ccc[1:25,], ccc2[46:50,])
    ccc.n2 = rbind(ccc2[1:45,], ccc[26:30,])
    xy = rbind(ccc.n1, ccc.n2,ccc.o1, ccc.o2)
  }
  
  if(n[1] == 20){ 
    cl = c(rep(1,15), rep(4, 5), 
           rep(2,55), rep(5, 5), 
           rep(3, 20))
    cl_new = c(rep(1,15), rep(-4, 5), 
               rep(2,55), rep(-5, 5), 
               rep(-3, 20))
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:20;  loca2 = 21:80;  loca3 = 81:100
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(20,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(60,0,sigma)
    mat[loca3,1] = simu_outlier(mat[loca3, 2], beta, distOut=2)
    mat = mat[,2:1] #change order
    
    # spatial coordinate
    ccc = mvrnorm(n[1], mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc2 = mvrnorm(n[2], mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
    ccc.o1 = mvrnorm(2, mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc.o2 = mvrnorm(18, mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2)) #18+2
    ccc.n1 = rbind(ccc[1:15,], ccc2[56:60,])
    ccc.n2 = rbind(ccc2[1:55,], ccc[16:20,])
    xy = rbind(ccc.n1, ccc.n2,ccc.o1, ccc.o2)
  }
  
  return(list(mat=mat, xy=xy, cl=cl, cl_new=cl_new, beta=beta))
}

simulateData_outlier_regOut <- function(n=c(40,40,20), beta=c(1.5,-1.5), sigma=0.05, 
                                     coordinate=c(1,1,-1,-1) 
)
{
  if(n[1] == 40){
    cl = c(rep(1,35), rep(4, 5), 
           rep(2,35), rep(5, 5), 
           rep(3, 20))
    cl_new = c(rep(1,35), rep(-4, 5), 
               rep(2,35), rep(-5, 5), 
               rep(-3, 20))
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:100
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
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
  }
  
  if(n[1] == 45){ 
    cl = c(rep(1,40), rep(4, 5), 
           rep(2,40), rep(5, 5), 
           rep(3, 10))
    cl_new = c(rep(1,40), rep(-4, 5), 
               rep(2,40), rep(-5, 5), 
               rep(-3, 10))
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:45;  loca2 = 46:90;  loca3 = 91:100
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(45,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(45,0,sigma)
    mat[loca3,1] = simu_outlier(mat[loca3, 2], beta, distOut=2)
    mat = mat[,2:1] #change order
    
    # spatial coordinate
    ccc = mvrnorm(n[1], mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc2 = mvrnorm(n[2], mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
    ccc.o1 = mvrnorm(n[3]/2, mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc.o2 = mvrnorm(n[3]/2, mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
    ccc.n1 = rbind(ccc[1:40,], ccc2[41:45,])
    ccc.n2 = rbind(ccc2[1:40,], ccc[41:45,])
    xy = rbind(ccc.n1, ccc.n2,ccc.o1, ccc.o2)
  }
  
  return(list(mat=mat, xy=xy, cl=cl, cl_new=cl_new, beta=beta))
}

simulateData_outlier_spaOut <- function(n=c(40,40,20), beta=c(1.5,-1.5), sigma=0.05, 
                                        coordinate=c(1,1,-1,-1),n_spa=10 
)
{
  if(n_spa == 10){
    cl = c(rep(1,35), rep(4, 5), 
           rep(2,35), rep(5, 5), 
           rep(3, 20))
    cl_new = c(rep(1,35), rep(-4, 5), 
               rep(2,35), rep(-5, 5), 
               rep(-3, 20))
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:100
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
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
  }
  
  if(n_spa == 20){ 
    cl = c(rep(1,30), rep(4, 10), 
           rep(2,30), rep(5, 10), 
           rep(3, 20))
    cl_new = c(rep(1,30), rep(-4, 10), 
               rep(2,30), rep(-5, 10), 
               rep(-3, 20))
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:100
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
    mat[loca3,1] = simu_outlier(mat[loca3, 2], beta, distOut=2)
    mat = mat[,2:1] #change order
    
    # spatial coordinate
    ccc = mvrnorm(n[1], mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc2 = mvrnorm(n[2], mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
    ccc.o1 = mvrnorm(n[3]/2, mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
    ccc.o2 = mvrnorm(n[3]/2, mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
    ccc.n1 = rbind(ccc[1:30,], ccc2[31:40,])
    ccc.n2 = rbind(ccc2[1:30,], ccc[31:40,])
    xy = rbind(ccc.n1, ccc.n2,ccc.o1, ccc.o2)
  }
  
  return(list(mat=mat, xy=xy, cl=cl, cl_new=cl_new, beta=beta))
}

# simulateData_outlier_shape <- function(sh='circle', n=c(40,40,20), beta=c(1.5,-1.5), sigma=0.05, 
#                                      coordinate=c(1,1,-1,-1) 
# )
# {
#   if(sh == 'circle'){
#     cl = c(rep(1,35), rep(4, 5), 
#            rep(2,35), rep(5, 5), 
#            rep(3, 20))

#     
#     mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
#     colnames(mat) = c('y','x')
#     loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:100
#     mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
#     mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
#     mat[loca3,1] = simu_outlier(mat[loca3, 2], beta, distOut=2)
#     mat = mat[,2:1] #change order
#     
#     # spatial coordinate
#     ccc = mvrnorm(n[1], mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
#     ccc2 = mvrnorm(n[2], mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
#     ccc.o1 = mvrnorm(n[3]/2, mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
#     ccc.o2 = mvrnorm(n[3]/2, mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.1, 2,2))
#     ccc.n1 = rbind(ccc[1:35,], ccc2[36:40,])
#     ccc.n2 = rbind(ccc2[1:35,], ccc[36:40,])
#     xy = rbind(ccc.n1, ccc.n2,ccc.o1, ccc.o2)
#   }
#   
#   ccc = mvrnorm(n[1], mu = c(coordinate[1],coordinate[2]), Sigma=diag(0.1, 2,2))
#   ccc2 = mvrnorm(n[2], mu = c(coordinate[3],coordinate[4]), Sigma=diag(0.5, 2,2))
#   xy = rbind(ccc,ccc2)
#   plot(xy)
#   
#   return(list(mat=mat, xy=xy, cl=cl, beta=beta))
# }


library('NonNorMvtDist')
simulateData_outlier_shape <- function(sh='gaus', n=c(40,40,20), beta=c(1.5,-1.5), sigma=0.05, 
                                       coordinate=c(1,1,-1,-1) 
)
{
  if(sh == 'gaus'){
    cl = c(rep(1,35), rep(4, 5), 
           rep(2,35), rep(5, 5), 
           rep(3, 20))
    cl_new = c(rep(1,35), rep(-4, 5), 
               rep(2,35), rep(-5, 5), 
               rep(-3, 20))
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:100
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
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
  }
  
  #------------------------
  # library('NonNorMvtDist') #example
  # tmp = rmvunif(n = 50, parm = 2, dim = 2)
  # tmp[,1] = tmp[,1] + 0.5
  # tmp[,2] = tmp[,2] + 0.5
  # plot(tmp)
  # 
  # tmp2 = rmvunif(n = 50, parm = 2, dim = 2)
  # tmp2[,1] = tmp2[,1] - 0.5
  # tmp2[,2] = tmp2[,2] - 0.5
  # plot(tmp2)
  # 
  # tmp_c = rbind(tmp,tmp2)
  # plot(tmp_c)
  mvruni1 <- function(n){
    tmp = rmvunif(n, parm = 2, dim = 2)
    tmp[,1] = tmp[,1] + 0.5
    tmp[,2] = tmp[,2] + 0.5
    return(tmp)
  }
  mvruni2 <- function(n){
    tmp = rmvunif(n, parm = 2, dim = 2)
    tmp[,1] = tmp[,1] - 0.5
    tmp[,2] = tmp[,2] - 0.5
    return(tmp)
  }
  
  #------------------------
  if(sh == 'unif'){
    cl = c(rep(1,35), rep(4, 5), 
           rep(2,35), rep(5, 5), 
           rep(3, 20))
    cl_new = c(rep(1,35), rep(-4, 5), 
               rep(2,35), rep(-5, 5), 
               rep(-3, 20))
    
    mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
    colnames(mat) = c('y','x')
    loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:100
    mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
    mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
    mat[loca3,1] = simu_outlier(mat[loca3, 2], beta, distOut=2)
    mat = mat[,2:1] #change order
    
    # spatial coordinate
    ccc = mvruni1(n=40)
    ccc2 = mvruni2(n=40)
    ccc.o1 = mvruni1(n=10)
    ccc.o2 = mvruni2(n=10)
    ccc.n1 = rbind(ccc[1:35,], ccc2[36:40,])
    ccc.n2 = rbind(ccc2[1:35,], ccc[36:40,])
    xy = rbind(ccc.n1, ccc.n2,ccc.o1, ccc.o2)
  }

  return(list(mat=mat, xy=xy, cl=cl, cl_new=cl_new, beta=beta))
}


simulateData_outlier_spaRadius <- function(n=c(40,40,20), beta=c(1.5,-1.5), sigma=0.05, 
                                     coordinate=c(1,1,-1,-1), rad=c(0.1,0.1)
)
{
  cl = c(rep(1,35), rep(4, 5), 
         rep(2,35), rep(5, 5), 
         rep(3, 20))
  cl_new = c(rep(1,35), rep(-4, 5), 
             rep(2,35), rep(-5, 5), 
             rep(-3, 20))
  
  mat = matrix(runif(sum(n)*2,-2,2), sum(n), 2)
  colnames(mat) = c('y','x')
  loca1 = 1:40;  loca2 = 41:80;  loca3 = 81:100
  mat[loca1,1] = beta[1] * mat[loca1,2] + rnorm(40,0,sigma)
  mat[loca2,1] = beta[2] * mat[loca2,2] + rnorm(40,0,sigma)
  mat[loca3,1] = simu_outlier(mat[loca3, 2], beta, distOut=2)
  mat = mat[,2:1] #change order
  
  # spatial coordinate
  ccc = mvrnorm(n[1], mu = c(coordinate[1],coordinate[2]), Sigma=diag(rad[1], 2,2))
  ccc2 = mvrnorm(n[2], mu = c(coordinate[3],coordinate[4]), Sigma=diag(rad[2], 2,2))
  ccc.o1 = mvrnorm(n[3]/2, mu = c(coordinate[1],coordinate[2]), Sigma=diag(rad[1], 2,2))
  ccc.o2 = mvrnorm(n[3]/2, mu = c(coordinate[3],coordinate[4]), Sigma=diag(rad[2], 2,2))
  ccc.n1 = rbind(ccc[1:35,], ccc2[36:40,])
  ccc.n2 = rbind(ccc2[1:35,], ccc[36:40,])
  xy = rbind(ccc.n1, ccc.n2,ccc.o1, ccc.o2)
  
  return(list(mat=mat, xy=xy, cl=cl, cl_new=cl_new, beta=beta))
}
