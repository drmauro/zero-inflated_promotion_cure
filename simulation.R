rm(list=ls(all=TRUE))
require(R2OpenBUGS) #library that links OPENBUGS to R
require(rootSolve)  #library for the multiroot command
require(VGAM)       #library to generate gamma-generalized values
require(coda)       #library to test the convergence of MCMCs
require(modeest)    #This package provides estimators of the mode
library(parallel)   #library to run in parallel


## ESTUDANDO O CENÁRIO
simFunc <- function(N){
  B<- 500    ### Numero de Estimativas
  set.seed(2016)
  ### Que serão utilizar no OPENBUGS
  cadeia<-800
  thining<-30
  Burnin<-300
  indice<-cadeia-Burnin
  
  
  #parameter setting
  beta1 = 0.5 
  beta2 = 0.5 
  beta3 = 1.5
  beta4 =   2
  beta5 =  -3
  beta6 =   1
  beta7 =  -2
  beta8 = 0.75
  varp=c(beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8)
  lime = 10
  
  #Declaring vectors and counters
  emcmc1 <-matrix(nrow=B,ncol=8)
  emcmc2 <-matrix(nrow=B,ncol=8)
  emcmc3 <-matrix(nrow=B,ncol=8)
  emle  <-matrix(nrow=B,ncol=8)
  old  =  options(digits=10)
  o<-1
  ite<-1
  pc1<-rep(0,times=B)
  pc2<-rep(0,times=B)
  pc3<-rep(0,times=B)
  pc4<-rep(0,times=B)
  pc5<-rep(0,times=B)
  pc6<-rep(0,times=B)
  pc7<-rep(0,times=B)
  pc8<-rep(0,times=B)
  pb1<-rep(0,times=B)
  pb2<-rep(0,times=B)
  pb3<-rep(0,times=B)
  pb4<-rep(0,times=B)
  pb5<-rep(0,times=B)
  pb6<-rep(0,times=B)
  pb7<-rep(0,times=B)
  pb8<-rep(0,times=B)
  censura = matrix(numeric(B),B,1)
  
  model <- function(N,varp){
    a       =  exp(varp[1])*exp(y*varp[2])
    lam     =  exp(varp[3])*exp(y*varp[4])
    
    ## PARAMETER WEIGHT IN ZEROS AND INFS: GAMMA ZERO AND GAMMA ONES WITH one (1) COVARIABLE
    gamm_zer = -log((exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8])))  
    tet      = -log((exp(varp[7])*exp(y*varp[8]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8])))
    
    U=runif(N)
    V=runif(N,exp(-gamm_zer),1-exp(-tet))
    y=numeric(N)
    for(i in 1:N){
      if(U[i]<=exp(-gamm_zer[i])){ 
        y[i] = 0
      }else 
        if(U[i]>exp(-gamm_zer[i]) & U[i] <= 1-exp(-tet[i])){
          y[i] = qweibull(-log(1-(V[i]-exp(-gamm_zer[i]))*(1-exp(-tet[i]))/(1-exp(-gamm_zer[i])-exp(-tet[i]) ) )/tet[i], shape=a[i], scale = lam[i], lower.tail = TRUE, log.p = FALSE)
        }else {      
          y[i]=+Inf
        }
    }
    return(y)
  }
  
  
  ## Declining the prior and posterior distributions in OPENBUGS ##
  modelo<-function () {
    
    for (i in 1:N)
    {
      f0[i] <- (alpha[i]/theta[i])*pow(t[i]/theta[i],alpha[i]-1)*exp(-pow(t[i]/theta[i],alpha[i]))
      F0[i] <- 1 - exp(-pow(t[i]/theta[i],alpha[i]))
      S0[i] <- exp(-pow(t[i]/theta[i],alpha[i]))
      
      f1[i] <- -(pow(p1[i],F0[i])/(1 - p1[i]))*log(p1[i])*f0[i]
      S1[i] <-  (pow(p1[i],F0[i]) - p1[i])/(1 - p1[i])
      
      L[i] <- ((1-step(t[i]-0.0000001))*p0[i]) +(step(t[i]-0.0000001)*pow((1-p0[i]-p1[i])*f1[i],d[i])*pow(p1[i]+(1-p0[i]-p1[i])*S1[i],1-d[i]))
      logL[i] <- log(L[i])
      
      alpha[i]<-exp(beta[1]+(x[i]*beta[2]))
      theta[i]<-exp(beta[3]+(x[i]*beta[4]))
      
      p0[i] <- exp(beta[5]+(x[i]*beta[6]))/(1+exp(beta[5]+(x[i]*beta[6]))+exp(beta[7]+(x[i]*beta[8])))
      p1[i] <- exp(beta[7]+(x[i]*beta[8]))/(1+exp(beta[5]+(x[i]*beta[6]))+exp(beta[7]+(x[i]*beta[8])))
      
      zeros[i] <- 0
      zeros[i]~dloglik(logL[i])
    }
    # Prior distributions
    for (j in 1:8) {
      beta[j]~dnorm(0,0.4)
    }
    
  }

  require(R2OpenBUGS) 
  require(rootSolve)  
  require(VGAM)       
  require(coda)       
  require(modeest) 
  
  GPE=function(varp)
  {
    a       =  exp(varp[1])*exp(y*varp[2])
    lam     =  exp(varp[3])*exp(y*varp[4])
    ## PARAMETER WEIGHT IN ZEROS AND INFS: GAMMA ZERO AND GAMMA ONES WITH A one(1) COVARIABLE
    gamm_zer = -log((exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8])))  
    tet      = -log((exp(varp[7])*exp(y*varp[8]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8])))
    
    f0       = exp(-gamm_zer)     
    fdw      = (exp(log(a)-log(lam)))*((exp(log(tempo)-log(lam)))^(a-1))*(exp(-(exp(log(tempo)-log(lam)))^a))
    Fdw      = 1-(exp(-(exp(log(tempo)-log(lam)))^a))
    
    fpop     = (1-exp(-gamm_zer)-exp(-tet))*(exp(-tet*Fdw)/(1-exp(-tet)))*tet*fdw
    Spop     = exp(-tet) + (1-exp(-gamm_zer)-exp(-tet))*((exp(-tet*Fdw)-exp(-tet))/(1-exp(-tet)))   
    f        = (fpop^delta)*(Spop^(1-delta))
    g        = ifelse(tempo==0, f0, f)
    adFunc   = sum(log(g))
    return(adFunc)
  } 
  
  while(o<=B)
  {
    #################-----------------------##################
    y=rbinom(N,1,0.5)
    tp = model(N,varp)
    tM = max(tp[tp<+Inf])
    Z=runif(N,0,tM)
    tempo=numeric(N)
    for (i in 1:N)
    {tempo[i] = min(tp[i],Z[i])
    }
    delta=numeric(N)
    for (i in 1:N) {
      if (tempo[i] < Z[i]){
        delta[i]=1
      }
    }
    
    t<-tempo
    d<-delta
    x<-y 
    
    fit=try(optim(varp,GPE,method="BFGS",hessian=TRUE,control=list(fnscale=-1)))
    estc = try(fit$par)
    Hc=try(fit$hessian); varic=try(-solve(Hc))
    
    if ( is.finite(estc[1]) & is.finite(estc[2]) & is.finite(estc[3]) & is.finite(estc[4]) & 
         is.finite(estc[5]) & is.finite(estc[6]) & is.finite(estc[7]) & is.finite(estc[8])) {
    if ( is.finite(varic[1,1]) & is.finite(varic[2,2]) & is.finite(varic[3,3]) & is.finite(varic[4,4]) & 
         is.finite(varic[5,5]) & is.finite(varic[6,6]) & is.finite(varic[7,7]) & is.finite(varic[8,8])){
    if ( varic[1,1] > 0 & varic[2,2] > 0 & varic[3,3] > 0 & varic[4,4] > 0 & 
         varic[5,5] > 0 & varic[6,6] > 0 & varic[7,7] > 0 & varic[8,8] > 0 &
         estc[1] < lime & estc[2] < lime & estc[3] < lime & estc[4] < lime &
         estc[5] < lime & estc[6] < lime & estc[7] < lime & estc[8] < lime ){
    
    #Declarando valores iniciais do OpenBUGS
    inits <- function(){
      list(beta = varp)
    }
    #################-----------------------##################
    
    data<-list("N", "t","d","x")
    #################-----------------------##################
    
    
    
    resultados<-try(bugs(data, inits, model.file = modelo,parameters = c("beta"),debug=FALSE ,n.burnin=Burnin,n.thin=thining,n.chains = 1, n.iter = cadeia,codaPkg = T))
    codaobject <- try(read.bugs(resultados))
    aux123<-try(codaobject[1,1][[1]])
    if(is.finite(aux123)){
    vb1<-NULL ; vb2<-NULL; vb3<-NULL; vb4<-NULL 
    vb5<-NULL; vb6<-NULL; vb7<-NULL; vb8<-NULL
    
    for (i in 1:indice) {
      vb1[i]<- codaobject[i,1][[1]]
      vb2[i]<- codaobject[i,2][[1]]
      vb3[i]<- codaobject[i,3][[1]]
      vb4[i]<- codaobject[i,4][[1]]
      vb5[i]<- codaobject[i,5][[1]]
      vb6[i]<- codaobject[i,6][[1]]
      vb7[i]<- codaobject[i,7][[1]]
      vb8[i]<- codaobject[i,8][[1]]
    }
    
    #par(mfrow=c(4,4))
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #ts.plot(vb1,ylab="B_1")
    #par(mai=c(0.77,0.5,0.1, 0.07))
    #acf(vb1,main="B_1")
    
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #ts.plot(vb2,ylab="B_2")
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #acf(vb2,main="B_2")
    
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #ts.plot(vb3,ylab="B_3")
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #acf(vb3,main="B_3")
    
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #ts.plot(vb4,ylab="B_4")
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #acf(vb4,main="B_4")
    
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #ts.plot(vb5,ylab="B_5")
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #acf(vb5,main="B_5")
    
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #ts.plot(vb6,ylab="B_6")
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #acf(vb6,main="B_6")
    
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #ts.plot(vb7,ylab="B_7")
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #acf(vb7,main="B_7")
    
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #ts.plot(vb8,ylab="B_8")
    #par(mai=c(0.5,0.5,0.1, 0.07))
    #acf(vb8,main="B_8")
    
    prb1i<-quantile(vb1, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
    prb1s<-quantile(vb1, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
    prb2i<-quantile(vb2, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
    prb2s<-quantile(vb2, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
    prb3i<-quantile(vb3, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
    prb3s<-quantile(vb3, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
    prb4i<-quantile(vb4, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
    prb4s<-quantile(vb4, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
    prb5i<-quantile(vb5, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
    prb5s<-quantile(vb5, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
    prb6i<-quantile(vb6, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
    prb6s<-quantile(vb6, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
    prb7i<-quantile(vb7, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
    prb7s<-quantile(vb7, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
    prb8i<-quantile(vb8, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
    prb8s<-quantile(vb8, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
 
    var.varp1 =varic[1,1]; Lvarp1 =estc[1]-1.96*sqrt(var.varp1); Uvarp1=estc[1]+1.96*sqrt(var.varp1)
    var.varp2 =varic[2,2]; Lvarp2 =estc[2]-1.96*sqrt(var.varp2); Uvarp2=estc[2]+1.96*sqrt(var.varp2)
    var.varp3 =varic[3,3]; Lvarp3 =estc[3]-1.96*sqrt(var.varp3); Uvarp3=estc[3]+1.96*sqrt(var.varp3)
    var.varp4 =varic[4,4]; Lvarp4 =estc[4]-1.96*sqrt(var.varp4); Uvarp4=estc[4]+1.96*sqrt(var.varp4)
    var.varp5 =varic[5,5]; Lvarp5 =estc[5]-1.96*sqrt(var.varp5); Uvarp5=estc[5]+1.96*sqrt(var.varp5)
    var.varp6 =varic[6,6]; Lvarp6 =estc[6]-1.96*sqrt(var.varp6); Uvarp6=estc[6]+1.96*sqrt(var.varp6)
    var.varp7 =varic[7,7]; Lvarp7 =estc[7]-1.96*sqrt(var.varp7); Uvarp7=estc[7]+1.96*sqrt(var.varp7)
    var.varp8 =varic[8,8]; Lvarp8 =estc[8]-1.96*sqrt(var.varp8); Uvarp8=estc[8]+1.96*sqrt(var.varp8)

    ge1<-abs(geweke.diag(vb1)$z[1])
    ge2<-abs(geweke.diag(vb2)$z[1])
    ge3<-abs(geweke.diag(vb3)$z[1])
    ge4<-abs(geweke.diag(vb4)$z[1])
    ge5<-abs(geweke.diag(vb5)$z[1])
    ge6<-abs(geweke.diag(vb6)$z[1])
    ge7<-abs(geweke.diag(vb7)$z[1])
    ge8<-abs(geweke.diag(vb8)$z[1])
    
    auxb1_1<- mlv(vb1, method = "grenander",p=4, all = TRUE, abc = FALSE)$M;  auxb1_2<- mlv(vb2, method = "grenander",p=4, all = TRUE, abc = FALSE)$M;  
    auxb1_3<- mlv(vb3, method = "grenander",p=4, all = TRUE, abc = FALSE)$M;  auxb1_4<- mlv(vb4, method = "grenander",p=4, all = TRUE, abc = FALSE)$M;  
    auxb1_5<- mlv(vb5, method = "grenander",p=4, all = TRUE, abc = FALSE)$M;  auxb1_6<- mlv(vb6, method = "grenander",p=4, all = TRUE, abc = FALSE)$M;  
    auxb1_7<- mlv(vb7, method = "grenander",p=4, all = TRUE, abc = FALSE)$M;  auxb1_8<- mlv(vb8, method = "grenander",p=4, all = TRUE, abc = FALSE)$M; 
    
    auxb2_1<- median(vb1);  auxb2_2<- median(vb2);  
    auxb2_3<- median(vb3);  auxb2_4<- median(vb4);  
    auxb2_5<- median(vb5);  auxb2_6<- median(vb6);  
    auxb2_7<- median(vb7);  auxb2_8<- median(vb8); 
  
    auxb3_1<- mean(vb1);  auxb3_2<- mean(vb2);  
    auxb3_3<- mean(vb3);  auxb3_4<- mean(vb4);  
    auxb3_5<- mean(vb5);  auxb3_6<- mean(vb6);  
    auxb3_7<- mean(vb7);  auxb3_8<- mean(vb8); 
    
    if (ge1<1.96  & ge2<1.96 & ge3<1.96  & ge4<1.96 & ge5<1.96& ge5<1.96  & ge6<1.96 & ge7<1.96) {  
    if (auxb1_1<lime & auxb1_2<lime & auxb1_3<lime & auxb1_4<lime & 
        auxb1_5<lime & auxb1_6<lime & auxb1_7<lime & auxb1_8<lime &
        auxb2_1<lime & auxb2_2<lime & auxb2_3<lime & auxb2_4<lime & 
        auxb2_5<lime & auxb2_6<lime & auxb2_7<lime & auxb2_8<lime &
        auxb3_1<lime & auxb3_2<lime & auxb3_3<lime & auxb3_4<lime & 
        auxb3_5<lime & auxb3_6<lime & auxb3_7<lime & auxb3_8<lime){   
    emcmc1[o,1] <- auxb1_1; emcmc1[o,2]<- auxb1_2;  emcmc1[o,3]<- auxb1_3; emcmc1[o,4] <- auxb1_4; 
    emcmc1[o,5] <- auxb1_5; emcmc1[o,6]<- auxb1_6;  emcmc1[o,7]<- auxb1_7; emcmc1[o,8] <- auxb1_8; 
    emcmc2[o,1] <- auxb2_1; emcmc2[o,2]<- auxb2_2;  emcmc2[o,3]<- auxb2_3; emcmc2[o,4] <- auxb2_4; 
    emcmc2[o,5] <- auxb2_5; emcmc2[o,6]<- auxb2_6;  emcmc2[o,7]<- auxb2_7; emcmc2[o,8] <- auxb2_8; 
    emcmc3[o,1] <- auxb3_1; emcmc3[o,2]<- auxb3_2;  emcmc3[o,3]<- auxb3_3; emcmc3[o,4] <- auxb3_4; 
    emcmc3[o,5] <- auxb3_5; emcmc3[o,6]<- auxb3_6;  emcmc3[o,7]<- auxb3_7; emcmc3[o,8] <- auxb3_8; 
    emle[o,1] <- estc[1]; emle[o,2] <- estc[2]; emle[o,3] <- estc[3]; emle[o,4] <- estc[4]; 
    emle[o,5] <- estc[5]; emle[o,6] <- estc[6]; emle[o,7] <- estc[7]; emle[o,8] <- estc[8]; 
    {
      if(prb1i<=beta1 & prb1s>= beta1) pb1[o]<-1
      if(prb2i<=beta2 & prb2s>= beta2) pb2[o]<-1
      if(prb3i<=beta3 & prb3s>= beta3) pb3[o]<-1
      if(prb4i<=beta4 & prb4s>= beta4) pb4[o]<-1
      if(prb5i<=beta5 & prb5s>= beta5) pb5[o]<-1
      if(prb6i<=beta6 & prb6s>= beta6) pb6[o]<-1
      if(prb7i<=beta7 & prb7s>= beta7) pb7[o]<-1
      if(prb8i<=beta8 & prb8s>= beta8) pb8[o]<-1
      if (beta1 >= Lvarp1 && beta1 <= Uvarp1) pc1[o]<-1
      if (beta2 >= Lvarp2 && beta2 <= Uvarp2) pc2[o]<-1
      if (beta3 >= Lvarp3 && beta3 <= Uvarp3) pc3[o]<-1
      if (beta4 >= Lvarp4 && beta4 <= Uvarp4) pc4[o]<-1
      if (beta5 >= Lvarp5 && beta5 <= Uvarp5) pc5[o]<-1
      if (beta6 >= Lvarp6 && beta6 <= Uvarp6) pc6[o]<-1
      if (beta7 >= Lvarp7 && beta7 <= Uvarp7) pc7[o]<-1
      if (beta8 >= Lvarp8 && beta8 <= Uvarp8) pc8[o]<-1
      
      censura[o,] = table(d)[1]/N 
    };
    cat(o,"     ",ite,"  ",round(sum(pc1)/o,3),"  ",round(sum(pc2)/o,3),"  ",round(sum(pc3)/o,3),"  ",round(sum(pc4)/o,3),"  ",round(sum(pc5)/o,3),"  ",round(sum(pc6)/o,3),"  ",round(sum(pc7)/o,3),"  ",round(sum(pc8)/o,3),"\n"); o<-(o+1);
    }
    }
    }}}}
    ite<-ite+1
  }
  
  
  ESTI1b  <-c(mean(emcmc1[,1]),mean(emcmc1[,2]),mean(emcmc1[,3]),mean(emcmc1[,4]),mean(emcmc1[,5]),mean(emcmc1[,6]),mean(emcmc1[,7]),mean(emcmc1[,8]))
  BIAS1b  <-c(mean(emcmc1[,1])-beta1, mean(emcmc1[,2])-beta2, mean(emcmc1[,3])-beta3, mean(emcmc1[,4])-beta4, mean(emcmc1[,5])-beta5, mean(emcmc1[,6])-beta6, mean(emcmc1[,7])-beta7, mean(emcmc1[,8])-beta8)
  RMSE1b  <-c(sqrt(mean((emcmc1[,1]-beta1)^2)),sqrt(mean((emcmc1[,2]-beta2)^2)),sqrt(mean((emcmc1[,3]-beta3)^2)),sqrt(mean((emcmc1[,4]-beta4)^2)),sqrt(mean((emcmc1[,5]-beta5)^2)),sqrt(mean((emcmc1[,6]-beta6)^2)),sqrt(mean((emcmc1[,7]-beta7)^2)),sqrt(mean((emcmc1[,8]-beta8)^2)))
  
  ESTI2b  <-c(mean(emcmc2[,1]),mean(emcmc2[,2]),mean(emcmc2[,3]),mean(emcmc2[,4]),mean(emcmc2[,5]),mean(emcmc2[,6]),mean(emcmc2[,7]),mean(emcmc2[,8]))
  BIAS2b  <-c(mean(emcmc2[,1])-beta1, mean(emcmc2[,2])-beta2, mean(emcmc2[,3])-beta3, mean(emcmc2[,4])-beta4, mean(emcmc2[,5])-beta5, mean(emcmc2[,6])-beta6, mean(emcmc2[,7])-beta7, mean(emcmc2[,8])-beta8)
  RMSE2b  <-c(sqrt(mean((emcmc2[,1]-beta1)^2)),sqrt(mean((emcmc2[,2]-beta2)^2)),sqrt(mean((emcmc2[,3]-beta3)^2)),sqrt(mean((emcmc2[,4]-beta4)^2)),sqrt(mean((emcmc2[,5]-beta5)^2)),sqrt(mean((emcmc2[,6]-beta6)^2)),sqrt(mean((emcmc2[,7]-beta7)^2)),sqrt(mean((emcmc2[,8]-beta8)^2)))
  
  ESTI3b  <-c(mean(emcmc3[,1]),mean(emcmc3[,2]),mean(emcmc3[,3]),mean(emcmc3[,4]),mean(emcmc3[,5]),mean(emcmc3[,6]),mean(emcmc3[,7]),mean(emcmc3[,8]))
  BIAS3b  <-c(mean(emcmc3[,1])-beta1, mean(emcmc3[,2])-beta2, mean(emcmc3[,3])-beta3, mean(emcmc3[,4])-beta4, mean(emcmc3[,5])-beta5, mean(emcmc3[,6])-beta6, mean(emcmc3[,7])-beta7, mean(emcmc3[,8])-beta8)
  RMSE3b  <-c(sqrt(mean((emcmc3[,1]-beta1)^2)),sqrt(mean((emcmc3[,2]-beta2)^2)),sqrt(mean((emcmc3[,3]-beta3)^2)),sqrt(mean((emcmc3[,4]-beta4)^2)),sqrt(mean((emcmc3[,5]-beta5)^2)),sqrt(mean((emcmc3[,6]-beta6)^2)),sqrt(mean((emcmc3[,7]-beta7)^2)),sqrt(mean((emcmc3[,8]-beta8)^2)))
  
  PCb    <-c(sum(pb1)/B,sum(pb2)/B,sum(pb3)/B,sum(pb4)/B,sum(pb5)/B,sum(pb6)/B,sum(pb7)/B,sum(pb8)/B)
  qtite  <-c(ite)
  s_cens = mean(censura[,1])

  ESTI1c <-c(mean(emle[,1]),mean(emle[,2]),mean(emle[,3]),mean(emle[,4]),mean(emle[,5]),mean(emle[,6]),mean(emle[,7]),mean(emle[,8]))
  BIAS1c <-c(mean(emle[,1])-beta1, mean(emle[,2])-beta2, mean(emle[,3])-beta3, mean(emle[,4])-beta4, mean(emle[,5])-beta5, mean(emle[,6])-beta6, mean(emle[,7])-beta7, mean(emle[,8])-beta8)
  RMSE1c <-c(sqrt(mean((emle[,1]-beta1)^2)),sqrt(mean((emle[,2]-beta2)^2)),sqrt(mean((emle[,3]-beta3)^2)),sqrt(mean((emle[,4]-beta4)^2)),sqrt(mean((emle[,5]-beta5)^2)),sqrt(mean((emle[,6]-beta6)^2)),sqrt(mean((emle[,7]-beta7)^2)),sqrt(mean((emle[,8]-beta8)^2)))
  PCc    <-c(sum(pc1)/B,sum(pc2)/B,sum(pc3)/B,sum(pc4)/B,sum(pc5)/B,sum(pc6)/B,sum(pc7)/B,sum(pc8)/B)

  return(list(ESTI1b  =  ESTI1b ,
              BIAS1b  =  BIAS1b ,
              RMSE1b  =  RMSE1b ,
              ESTI2b  =  ESTI2b ,
              BIAS2b  =  BIAS2b ,
              RMSE2b  =  RMSE2b ,
              ESTI3b  =  ESTI3b ,
              BIAS3b  =  BIAS3b ,
              RMSE3b  =  RMSE3b ,
              PCb     =  PCb    , 
              qtite   =  qtite  ,
              s_cens  =  s_cens ,
              ESTI1c  =  ESTI1c ,
              BIAS1c  =  BIAS1c ,
              RMSE1c  =  RMSE1c ,
              PCc     =  PCc   ) )
}


numCores <- 5
cl <- makeCluster(numCores)
finalResults1 <- parLapply(cl, c(100,250,500,750,1000), simFunc)
stopCluster(cl)

save.image("modelo2c1.RData") 