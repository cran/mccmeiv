.m4noa1=function(capm, y,sv,w,xs,z,len.xs,len.sv,len.z, no.par.one,no.par.two,startingvalues,naive.estimate.beta){
  ##
  ## Our proposed method 
  ###### estimation of the first likelihood
  #
  #out300=glm(w~tstar+sv+z, family=binomial)
  
  nlglikm4=function(par){
    ####Extract parameters from par object
    par.int=par[1]
    par.xs=par[2:(1+len.xs)]
    par.sv=par[(2+len.xs):(1+len.xs+len.sv)]
    par.z=par[(2+len.xs+len.sv):(1+len.xs+len.sv+len.z)]
    par.eta=par[(no.par.one)]
    #
    par.beta=par[(no.par.one+1):(no.par.one+no.par.two)]
    beta.g=par.beta[1]
    beta.z=par.beta[2:no.par.two]
    par.set=c(par.int,par.xs,par.sv,par.z,par.eta)
    #
    ## prpi=pr(X=1|X^*, Z) # a logistic model is assumed
    sum.xs=as.vector(as.matrix(xs)%*%par.xs)
    sum.sv=as.vector(as.matrix(sv)%*%par.sv)
    sum.z=as.vector(as.matrix(z)%*%par.z)
    prt1.0=1/(1+exp(-par.int-sum.xs-sum.sv-sum.z))
    prt1.1=1/(1+exp(-par.int-beta.g-sum.xs-sum.sv-sum.z))
    
    #In the event that we get an extremely high or extremely low
    #probability of treatment, set it to a specific value
    prt1.0[prt1.0>0.99999]=0.99999
    prt1.0[prt1.0<0.00001]=0.00001
    prt1.1[prt1.1>0.99999]=0.99999
    prt1.1[prt1.1<0.00001]=0.00001
    #
    ## alpha.0.1=pr(W=1|X=0, Z)     
    alpha.0.1=1/(1+exp(-par.eta[1]) )
    ##   alpha.1.1=pr(W=0|X=1, Z)=0
    
    ##
    # prw1.0=pr(W=1|T^*, Z,Y=0)
    # prw1.1=pr(W=1|T^*, Z,Y=1)
    prw1.0=alpha.0.1+(1-alpha.0.1)*prt1.0
    prw1.1=alpha.0.1+(1-alpha.0.1)*prt1.1
    
    prw1.0[prw1.0>0.99999]=0.99999
    prw1.0[prw1.0<0.00001]=0.00001
    prw1.1[prw1.1>0.99999]=0.99999
    prw1.1[prw1.1<0.00001]=0.00001
    
    
  
    ########################
    #This section is used for the beta formulas
    #g2=log(1-prt1.0+exp(par[7])*prt1.0)
    
    exp.g2=(1-prt1.0+exp(beta.g)*prt1.0)
    
    
    ################ conditional probability calculation 
    sum.xs=as.vector(as.matrix(z)%*%beta.z)
    term1.2=exp.g2*exp(sum.xs) 
    term2.2=matrix(term1.2, ncol=(capm+1), byrow=T)
    term3.2=rep(apply(term2.2, 1, sum), each=(capm+1))
    cond.prob.2=term1.2/term3.2
    
    cond.prob.2[cond.prob.2<0.00001]=0.00001
    #########################
    
    neglk=- sum (y*log(cond.prob.2)+
                   (y)*(log(prw1.1)*w+(1-w)*log(1-prw1.1))+
                   (1-y)*(log(prw1.0)*w+(1-w)*log(1-prw1.0)))
    
    return(neglk)}
  
  start=c(startingvalues,naive.estimate.beta)
  out2=optim(start, nlglikm4, method="L-BFGS-B", control=list(maxit=1500), hessian=T)
  par.set.m4=out2$par[1:no.par.one]
  
  beta1=no.par.one+1
  betaend=no.par.one+no.par.two
  beta=out2$par[beta1:betaend]
  
  hmat=out2$hessian
  var.cov=solve(hmat)
  beta.sd=sqrt(diag(var.cov)[beta1:betaend])
  par.set.sd=sqrt(diag(var.cov))
  
  ##
  ########################################################################### 
  
  
  beta.m4=beta
  var.cov.m4=var.cov
  beta.sd.m4=beta.sd
  
  result=list(par.set.m4=par.set.m4,par.set.sd.m4=par.set.sd,beta.m4=beta.m4,
              var.cov.m4=var.cov.m4,beta.sd.m4=beta.sd.m4)
  return(result)
  }