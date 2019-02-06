#This calculates the score functions "u" used in estimating the middle term for the standard
#error calculations

#Perhaps use data from nlglik2.r rather than recalculating data
#
.part3noa0=function(beta,par,y,sv,w,xs,z,prt1_1,alpha.1.1,prw1_1,len.xs,len.sv,len.z,no.par.one,no.par.two,n,capm,nallpara){
  

  ##Next we calculate the standard errors of all of the parameters
  ##First we need to calculate the "middle term" of the asymptotic variance
  ##It is composed of the U equations
  ##We need the output from the first analysis
  allpara=c(par, beta)
  no.par.one=length(par) 
  no.par.two=length(beta)
  nallpara=length(allpara) #This is used to determine the size of the uvec 
  
  beta.g=beta[1]
  beta.z=beta[2:no.par.two]
  
  exp.g.w1.1= (exp(beta.g)*(1-alpha.1.1)*prt1_1)/prw1_1
  exp.g.w0.1= (exp(beta.g)*alpha.1.1*prt1_1+(1-prt1_1))/(1-prw1_1)
  
  deriv.beta1.g.w1.1=1
  
  deriv.beta1.g.w0.1=(exp(beta.g)*alpha.1.1*prt1_1)/(exp(beta.g)*alpha.1.1*prt1_1+(1-prt1_1))
  
  deriv.beta1.g.1=deriv.beta1.g.w1.1*w+(1-w)*deriv.beta1.g.w0.1
  
  exp.g.1=exp.g.w1.1*w+(1-w)*exp.g.w0.1
  
  ################ conditional probability calculation
  sum.xs=as.vector(as.matrix(z)%*%beta.z)
  term1.2=exp.g.1*exp(sum.xs)
  term2.2=matrix(term1.2, ncol=(capm+1), byrow=T)
  term3.2=rep(apply(term2.2, 1, sum), each=(capm+1))
  cond.prob.2=term1.2/term3.2
  
  #
  mid.term.1=0
  
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
  
  ##After calculating the uvec, we can now calculate the middle term 
  #for asymptotic variance
  ##mid.term
  
  mid.term1.1=t(uvec.1)%*%uvec.1
  mid.term1.1=mid.term1.1/n
  mid.term.1=mid.term1.1
  
  return(list(uvec.1=uvec.1,mid.term.1=mid.term.1))
}