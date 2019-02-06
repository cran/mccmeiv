####### Derivative calculation using numerical approach #####
##After calculating the middle term, we now need to calculate A hat
##Once we have A hat we then  now calculate asymptotic variance
################
.part4noa0=function(mid.term.1,y,sv,w,xs,z,allpara,len.xs,len.sv,len.z,no.par.one,no.par.two,n,capm,nallpara){

  f2=function(allpara){
  #Convert data into vectors
   
  #extract parameters
  par=allpara[1:no.par.one]
  beta= allpara[(no.par.one+1):nallpara]
  
  #extract parameters
  par.int=par[1]
  par.xs=par[2:(1+len.xs)]
  par.sv=par[(2+len.xs):(1+len.xs+len.sv)]
  par.z=par[(2+len.xs+len.sv):(1+len.xs+len.sv+len.z)]
  par.eta=par[no.par.one]

  #
  ## prpi=pr(X=1|X^*, Z) # a logistic model is assumed 
  sum.xs=as.vector(as.matrix(xs)%*%par.xs)
  sum.sv=as.vector(as.matrix(sv)%*%par.sv)
  sum.z=as.vector(as.matrix(z)%*%par.z)
  prt1_1=1/(1+exp(-par.int-sum.xs-sum.sv-sum.z))
  #
  ## alpha.0.1=pr(W=1|X=0, Z)=0     
    ##   alpha.1.1=pr(W=0|X=1, Z)
  alpha.1.1=1/(1+exp(-par.eta[1]) )
  
  #
  # prw.1=pr(W=1|X^*, Z)
  prw1_1=(1-alpha.1.1)*prt1_1
  
  #Function used to find optimal values of beta


    beta.g=beta[1]
    beta.z=beta[2:no.par.two]

    #calculate terms for derivatives
    exp.g.w1.1= (exp(beta.g)*(1-alpha.1.1)*prt1_1)/prw1_1
    exp.g.w0.1= (exp(beta.g)*alpha.1.1*prt1_1+(1-prt1_1))/(1-prw1_1)
    
    deriv.beta1.g.w1.1=1
    
    deriv.beta1.g.w0.1=(exp(beta.g)*alpha.1.1*prt1_1)/(exp(beta.g)*alpha.1.1*prt1_1+(1-0)*(1-prt1_1))
    
    deriv.beta1.g.1=deriv.beta1.g.w1.1*w+(1-w)*deriv.beta1.g.w0.1
    
    exp.g.1=exp.g.w1.1*w+(1-w)*exp.g.w0.1
    
    ################ conditional probability calculation
    sum.xs=as.vector(as.matrix(z)%*%beta.z)
    term1.2=exp.g.1*exp(sum.xs) 
    term2.2=matrix(term1.2, ncol=(capm+1), byrow=T)
    term3.2=rep(apply(term2.2, 1, sum), each=(capm+1))
    cond.prob.2=term1.2/term3.2
    
    ##Here the estimates for the U equations are estimates
    uvec.1=matrix(NA,nrow=n,ncol=nallpara,byrow=T)
    dim(uvec.1)
    
    #U equations for gamma
    #Need terms involving the instrument and the stratification variables
    #Will need for loops because we don't know the number of strat variables and instruments
    term11=(
      (1-y)*(1-alpha.1.1)*prt1_1*
        (1-prt1_1)*(w-prw1_1)/(prw1_1*(1-prw1_1))
    )
    
    term1.mat=matrix(term11, ncol=(capm+1), byrow=T)
    dim(term1.mat)
    summary(term1.mat)
    #instruments
    uvec.1[,1]=apply(term1.mat,1,sum)
    for (i in 1:len.xs){ 
      uvecxs=term11*xs[,i]
      uvecxs.mat=matrix(uvecxs,ncol=(capm+1),byrow=T)
      uvec.1[,(i+1)]=apply(uvecxs.mat,1,sum)
    }
    
    #stratification variables
    for (i in 1:len.sv){ 
      uvecsv=term11*sv[,i]
      uvecsv.mat=matrix(uvecsv,ncol=(capm+1),byrow=T)
      uvec.1[,(1+len.xs+i)]=apply(uvecsv.mat,1,sum)
    }
    
    #prognostic factors
    for (i in 1:len.z){ 
      uvecz=term11*z[,i]
      uvecz.mat=matrix(uvecz,ncol=(capm+1),byrow=T)
      uvec.1[,(1+len.xs+len.sv+i)]=apply(uvecz.mat,1,sum)
    }
    
    #U equations for eta
    term2= (1-y)*(w-prw1_1)/(prw1_1*(1-prw1_1))
    
    uveceta2= term2*( -alpha.1.1*(1-alpha.1.1)*prt1_1) 
    uveceta2.mat=matrix(uveceta2,ncol=(capm+1),byrow=T)
    uvec.1[,(1+len.xs+len.sv+len.z+1)]=apply(uveceta2.mat,1,sum)
    
    #U equations for beta
    newterm1.1=y-cond.prob.2
    
    mult.1=(newterm1.1)*deriv.beta1.g.1
    mult.mat.1=matrix(mult.1,ncol=(capm+1),byrow=T)
    uvec.1[,(1+len.xs+len.sv+len.z+2)]=apply(mult.mat.1,1,sum)
    
    for (i in 1:len.z){
      uvecbetaz=newterm1.1*z[,i]
      uvecbetaz.mat=matrix(uvecbetaz,ncol=(capm+1),byrow=T)
      uvec.1[,(1+len.xs+len.sv+len.z+2+i)]=apply(uvecbetaz.mat,1,sum)
    }
    
    return(apply(uvec.1, 2, sum))
    
    #End function
  }
 
  
  #
  #
  #
  ##Finally calculate the A hat, then calculate the asymptotic variance
  deriv.inv=ginv(jacobian(f2, allpara)/n)
  deriv.1=deriv.inv
  var.cov.1=deriv.inv%*%mid.term.1%*%t(deriv.inv)/n
  theta.sd.1=sqrt(diag(var.cov.1))
   
##########
##########


return(list(theta.sd.1=theta.sd.1,var.cov.1=var.cov.1,deriv.1=deriv.1))

#End function
  }
