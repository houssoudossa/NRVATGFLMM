NRVAT <-
function(n, gamm, heritability, a1=1, b1=25, Ped, covariates, geno, kin2){
  
  
  ##############################################
  # Function for getting the Weigth of the SNP #
  ##############################################
  Generate.W=function(a1,b1,p)
  {
    N=length(p)
    return(diag(dbeta(p,a1,b1),N,N))
  }
  #diag(dbeta(maf, w_a, w_b))
  
  Estim.Gamm = function(X,Y)
  {
    X1=X[,2]
    X2=X[,3]
    fit = glm(Y~X1+X2,family=binomial)
    gamm.hat=summary(fit)$coef[,1]
    return(gamm.hat)
  }
  
  ###################################################
  ##################### Estime bien h^2 #############
  ###################################################
  like_general<-function(heritability,kin2,mu.hat,Y){
    I=I
    n=n
    
    #X = Generate.X(n)
    #mu = Compute.mu(X,gamm)
    #Y = Generate.Y(I,n,mu,heritability,kin2)
    #gamm.hat = Estim.Gamm(X,Y)
    #mu.hat = Compute.mu(X,gamm.hat)
    
    llikelihood = rep(0,I)
    for(i in 1:I)
    {
      indices.fam=split(1:sum(n), rep(1:I, n))
      Y.fam = Y[indices.fam[[i]]]
      mu.fam = mu.hat[indices.fam[[i]]]
      
      Sigma=heritability*2*kin2[indices.fam[[i]],indices.fam[[i]]]+(1-heritability)*diag(rep(1,n[i]))
      
      for(j in 1:(n[i]-1))
      {
        for(k in (j+1):n[i])
        {
          ### Using MVTNORM ######
          ########################
          #Sigma1=Sigma[c(j,k),c(j,k)]
          
          #copCDF.u1.u2<- pmvnorm(upper=c(qnorm(mu.fam[j]),qnorm(mu.fam[k])),sigma=Sigma[c(j,k),c(j,k)])
          
          ### Using Copula #######
          ########################
          
          cop <- BiCop(family = 1, Sigma[j,k])
          copCDF.u1.u2<- BiCopCDF(mu.fam[j],mu.fam[k],cop)
          
          if ((Y.fam[j]==1)&&(Y.fam[k]==1)) lik.fam = log(copCDF.u1.u2) 
          if ((Y.fam[j]==0)&&(Y.fam[k]==1)) lik.fam = log(mu.fam[k]-copCDF.u1.u2)
          if ((Y.fam[j]==1)&&(Y.fam[k]==0)) lik.fam = log(mu.fam[j]-copCDF.u1.u2)
          if ((Y.fam[j]==0)&&(Y.fam[k]==0)) lik.fam = log(1-mu.fam[j]-mu.fam[k]+copCDF.u1.u2) 
          
          lik.fam=sum(lik.fam)
          
        }
        
      }
      
      llikelihood[[i]] = lik.fam
    }
    #return(sum(llikelihood))
    return(-sum(llikelihood))
  }
  
  ############################
  ## Date : 26 - 04 - 2019  ##
  ## Auteur : Roland        ##
  ############################
  
  #############################################
  ## Computation of Mu under Null hypothesis ##
  #############################################
  
  Compute.mu = function(X,gamm) {as.vector((exp(as.matrix(X)%*%gamm)/(1+exp(as.matrix(X)%*%gamm))))}
  
  ####################################################
  ## Computation of Sigma with the Heritability h^2 ##
  ####################################################
  
  compute.Sigma = function(heritability,i) {heritability*rep(1,n[i])%*%t(rep(1,n[i]))+(1-heritability)*diag(rep(1,n[i]))}
  
  ################################
  ## Computation of The matrice ## 
  ################################
  
  Compute.V = function(mu,I,n,heritability,kin2)
    
  {
    V = NULL
    for(i in 1:I){
      indices.fam=split(1:sum(n), rep(1:I, n))
      
      mu.fam = mu[indices.fam[[i]]]
      V.fam = diag(mu.fam*(1-mu.fam))
      
      Sigma=heritability*2*kin2[indices.fam[[i]],indices.fam[[i]]]+(1-heritability)*diag(rep(1,n[i]))
      
      for(j in 1:(n[i]-1))
      {
        for(k in (j+1):n[i])
        {
          ############################
          ### Using MVTNORM ##########
          ############################
          #Sigma.2=Sigma[c(j,k),c(j,k)]
          #V.fam[j,k] = pmvnorm(upper=c(qnorm(mu.fam[j]),qnorm(mu.fam[k])),sigma=Sigma[c(j,k),c(j,k)])- mu.fam[j]*mu.fam[k]
          
          #V.fam[j,k]=pmvnorm(upper=c(qnorm(mu.fam[j]),qnorm(mu.fam[k])),corr=Sigma.2) - mu.fam[j]*mu.fam[k]
          
          ############################
          ######## Using Copula ######
          ############################
          
          cop <- BiCop(family = 1, Sigma[j,k])
          
          copCDF.u1.u2<- BiCopCDF(mu.fam[j],mu.fam[k],cop)
          
          V.fam[j,k]<- copCDF.u1.u2-mu.fam[j]*mu.fam[k]
          
          V.fam[k,j] <- V.fam[j,k]
        }
      }
      V[[i]] = V.fam
    }
    bdiag(V)
  }
  
  ###########################################
  ## Computation of A (Fisher Information) ##
  ###########################################
  
  Compute.A = function(X,I,mu) {(t(as.matrix(X))%*%diag(mu*(1-mu))%*%as.matrix(X))/(I)}
  
  #######################
  ## Computation of B ###
  #######################
  
  Compute.B = function(X,V,I) {(t(as.matrix(X))%*%V%*%as.matrix(X))/(I)}
  
  #########################################
  ## Computation of D (derivative of Mu) ##
  #########################################
  
  Compute.D = function(mu,X,n){
    D = NULL
    I=length(n)
    for(i in 1:I)
    {
      indices.fam=split(1:sum(n), rep(1:I, n))
      
      mu.fam = mu[indices.fam[[i]]]
      X.fam = as.matrix(X)[indices.fam[[i]],]
      result = diag(mu.fam*(1-mu.fam))%*%X.fam
      D[[i]] = result
    }
    #As D is liste, we transforme it as matrice
    D=do.call(rbind,D)
    D
    #bdiag(D)
  }
  
  
  
  #####################
  ## Kernel fonctions##
  #####################
  
  K1_Help= function(x,y){
    # Helper function for 2 way interaction kernel
    p = length(x)
    a = x*y
    b = cumsum(a)
    return(sum(a[-1]*b[-p]))
  }
  
  
  call_Kernel_IBS<-function(Z,n,p){
    
    #Kernel_IBS(double * Z, int * pn, int * pp, double * Kernel)
    K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
    temp<-.C("Kernel_IBS",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.double(as.vector(K)))[[4]]
    matrix(temp,nrow=n)
  }
  
  
  
  call_Kernel_IBS_Weight<-function(Z,n,p,weights){
    
    #Kernel_IBS_Weight(int * Z, int * pn, int * pp, int *UseGivenWeight ,  double * weight, double * Kernel)
    given_weight = 1;
    if( is.null(weights)){
      weights = rep(0,p);
      given_weight = 0;
    } 
    else {
      # change!!
      weights<-weights^2;
    }
    K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
    temp<-.C("Kernel_IBS_Weight",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.integer(given_weight),
             as.double(weights),as.double(as.vector(K)))[[6]]
    matrix(temp,nrow=n)
  }
  
  
  call_Kernel_2wayIX<-function(Z,n,p){
    
    #Kernel_IBS(double * Z, int * pn, int * pp, double * Kernel)
    K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
    temp<-.C("Kernel_2wayIX",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.double(as.vector(K)))[[4]]
    matrix(temp,nrow=n)
  }
  
  ################################################################################
  # GAUSSIAN KERNEL rho > 0 
  ################################################################################
  
  kernel.gaussian <- function(x, rho)
  {
    out <- exp(-1*as.matrix(dist(x) ^ 2)/rho)
    return(out)
  }
  ################################################################################
  # POLYNOMIAL KERNEL
  ################################################################################
  
  kernel.polynomial <- function(x, rho, gamma, d)
  {
    return((rho * x %*% t(x) + gamma) ^ d)
  }
  
  ################################################################################
  # Integrating the differents KernelL functions
  ################################################################################
  KNL = function(Z, kernel, weights,n,m,rho,gamma,d){
    ## Add linear, linear.W and quadratic.w
    if (kernel == "linear") {
      K = Z%*%t(Z)
    }
    
    if (kernel == "linear.w") {
      K = Z%*%weights%*%t(Z)
    }
    
    if (kernel == "gauss") {
      K = exp(-1*as.matrix(dist(Z) ^ 2)/rho)
    }
    
    if (kernel == "gauss.w") {
      K = exp(-1*as.matrix(dist(Z %*%weights%*%t(Z))^ 2)/rho)
    }
    
    if (kernel == "poly") {
      K = (rho * Z %*% t(Z) + gamma)^ d
    }
    
    if (kernel == "poly.w") {
      K = (rho * Z %*%weights%*%t(Z) + gamma)^ d
    }
    
    if (kernel == "quadratic") {
      K = (Z%*%t(Z)+1)**2
    }
    
    if (kernel == "quadratic.w") {
      K = (Z%*%weights%*%t(Z)+1)**2
    }
    
    if (kernel == "IBS") {
      K = call_Kernel_IBS(Z,n,m)
    }
    if (kernel == "IBS.weighted") {
      
      K = call_Kernel_IBS_Weight(Z,n,m,weights)
    }
    if (kernel == "2wayIX") {
      K = call_Kernel_2wayIX(Z,n,m)
    }  
    if (kernel == "IBS.weighted_OLD") {
      #K = matrix(nrow = n, ncol = n)
      if (is.null(weights)) {
        qs = apply(Z, 2, mean)/(2)
        weights = 1/sqrt(qs)
      } else {
        weights<-weights^2
      }
      K1 = matrix(nrow =n, ncol = n)
      for (i in 1:n) {
        K1[i,] = apply(abs(t(Z)-Z[i,])*weights,2, sum)
      }
      K= 1-(K1)/(2*sum(weights))
    }
    
    if (kernel == "IBS_OLD") {
      K1=matrix(nrow=n,ncol=n)
      for (i in 1:n) {
        K1[i,] = apply(abs(t(Z)-Z[i,]),2, sum)
      }
      K = (2*m-K1)/(2*m)
    }
    if (kernel == "2wayIX_OLD") {
      K = 1+Z%*%t(Z)
      N1=  matrix(nrow = n, ncol = n)
      for (i in 1:n){
        for (j in i:n){
          N1[j,i] = N1[i,j] = K1_Help(Z[i,], Z[j,])
        }
      }
      K = K+N1
    }
    return(K)
    
  }
  
  #############################################################
  ######### Call The existing packages using library ##########
  #############################################################
  library(mvtnorm)
  library(CompQuadForm)
  library(VineCopula)
  library(copula)
  library(SKAT)
  library(Matrix)
  library(MASS)
  
  #n=rep(c(3,4,8),40)# n (number of individuals in each family) is different from one family to another.
  I=length(n) # Number of families
  
  #################### Parameters ##############
  #############################################
  X = covariates[,-c(1,2)]
  mu = Compute.mu(X,gamm)
  Y = Ped$trait
  gamm.hat = Estim.Gamm(X,Y)
  mu.hat = Compute.mu(X,gamm.hat)
  geno<-geno[, -c(1,2)]
  ############################## 
  #Getting Weight for the SNP 
  p=apply(geno, 2, function(x) {sum(x, na.rm = T)/((length(x) - sum(is.na(x))) * 2)}) #p1=runif(ncol(geno),a,b)
  
  W=Generate.W(a1,b1,p) #W=Generate.W(a1,b1,p1)
  ############################## 
  ###############################################
  ## Estimate of the dependance parameter (h^2)##
  ###############################################
  heritability_est=optim(par = heritability, fn = like_general, kin2=kin2, mu.hat=mu.hat,Y=Y, method = "Brent", lower = heritability, upper = heritability + 0.01)$par
  
  V.hat = Compute.V(mu.hat,I,n,heritability_est,kin2)
  #V.hat = Compute.V(mu.hat,I,n,heritability,kin2)
  D = Compute.D(mu.hat,X,n)
  A=Compute.A(X,I,mu.hat) 
  B=Compute.B(X,V.hat,I)
  I.A=solve(A)
  M1=V.hat%*%as.matrix(X)%*%I.A%*%t(D)
  M2= D%*%solve(A)%*%B%*%I.A%*%t(D)
  
  ### Computation of the Corrected Matrix
  
  MAtt.Cov= V.hat - (1/I)*t(M1)-(1/I)*M1+(1/I)*M2
  
  ### Eigen Decomposition of the Corrected Matrix
  sigmaY = eigen(MAtt.Cov)
  U = sigmaY$vectors
  delta.1 = sigmaY$values
  delta.1[which(delta.1<0)]=0
  sqrt.delta.1 = sqrt(delta.1)
  sqrt.sigmaY = U%*%diag(sqrt.delta.1,sum(n),sum(n))%*%t(U)
  
  ### Computation of the Linear Kernel Functions
  AL=KNL(Z=as.matrix(geno), kernel="linear.w",weights=W, n=sum(n),m=ncol(as.matrix(geno)),rho=ncol(as.matrix(geno)),gamma=1,d=2)
  KL=sqrt.sigmaY%*%AL%*%sqrt.sigmaY 
  eig=eigen(KL)
  DL=eig$values
  
  ### Computation of the P.Values
  S_obsL=t(Y-mu.hat) %*% AL %*%(Y-mu.hat)
  p.valueL=davies(S_obsL, DL)$Qq
  
  
  ### Computation of the Quadratic Kernel Functions
  AQ=KNL(Z=as.matrix(geno), kernel="quadratic.w",weights=W, n=sum(n),m=ncol(as.matrix(geno)),rho=ncol(as.matrix(geno)),gamma=1,d=2)
  KQ=sqrt.sigmaY%*%AQ%*%sqrt.sigmaY 
  eig=eigen(KQ)
  DQ=eig$values
  
  ### Computation of the P.Values
  S_obsQ=t(Y-mu.hat) %*% AQ %*%(Y-mu.hat)
  p.valueQ=davies(S_obsQ, DQ)$Qq
  
  
  ### Computation of the IBS Kernel Functions
  AIB=KNL(Z=as.matrix(geno), kernel="IBS.weighted",weights=W, n=sum(n),m=ncol(as.matrix(geno)),rho=ncol(as.matrix(geno)),gamma=1,d=2)
  KIB=sqrt.sigmaY%*%AIB%*%sqrt.sigmaY 
  eig=eigen(KIB)
  DIB=eig$values
  #DQ[which(DQ<0)]=0
  ### Computation of the P.Values
  S_obsIB=t(Y-mu.hat) %*% AIB %*%(Y-mu.hat)
  ###Get the Quadratic Statistic and P.value
  p.valueIB=davies(S_obsIB, DIB)$Qq
  
  
  ### Computation of the Gaussien Kernel Functions
  AG=KNL(Z=as.matrix(geno), kernel="gauss.w",weights=W, n=sum(n),m=ncol(as.matrix(geno)),rho=ncol(as.matrix(geno)),gamma=1,d=2)
  KG=sqrt.sigmaY%*%AG%*%sqrt.sigmaY 
  eig=eigen(KG)
  DG=eig$values
  #DG[which(DG<0)]=0
  ### Computation of the P.Values
  S_obsG=t(Y-mu.hat) %*% AG %*%(Y-mu.hat)
  ###Get the Gaussien Statistic and P.values
  p.valueG=davies(S_obsG, DG)$Qq
  
  
  ### Computation of the Polynomial Kernel Functions
  AP=KNL(Z=as.matrix(geno), kernel="poly.w",weights=W, n=sum(n),m=ncol(as.matrix(geno)),rho=ncol(as.matrix(geno)),gamma=1,d=2)
  KP=sqrt.sigmaY%*%AP%*%sqrt.sigmaY 
  eig=eigen(KP)
  DP=eig$values
  #DG[which(DG<0)]=0
  ### Computation of the P.Values
  S_obsP=t(Y-mu.hat) %*% AP %*%(Y-mu.hat)
  ###Get the Gaussien Statistic and P.values
  p.valueP=davies(S_obsP, DP)$Qq
  
  #### Table of Results######
  
  NRVAT=cbind(gamm.hat[1], gamm.hat[2],gamm.hat[3], heritability_est, S_obsL, p.valueL, S_obsQ, p.valueQ, S_obsIB, p.valueIB, S_obsG, p.valueG, S_obsP, p.valueP)
  colnames(NRVAT)=c("Intercept", "gam1", "gam2", "h.S_estim", "NRVAT_S_obsL", "NRVAT_p.valueL", "NRVAT_S_obsQ", "NRVAT_p.valueQ", "NRVAT_S_obsIB", "NRVAT_p.valueIB", "NRVAT_S_obsG", "NRVAT_p.valueG", "NRVAT_S_obsP", "NRVAT_p.valueP")
  rownames(NRVAT) = NULL
  return(as.data.frame(NRVAT))
}
