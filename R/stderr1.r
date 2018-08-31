#This calculates the score functions "u" used in estimating the middle term for the standard
#error calculations

#Perhaps use data from nlglik2.r rather than recalculating data
#
.part3=function(beta,par,y,sv,w,xs,z,prt1_1,alpha.0.1,alpha.1.1,prw1_1,len.xs,len.sv,len.z,no.par.one,no.par.two,n,capm,nallpara){
  

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
  
  exp.g.w1.1= (exp(beta.g)*(1-alpha.1.1)*prt1_1 +alpha.0.1*(1-prt1_1))/ prw1_1
  exp.g.w0.1= (exp(beta.g)*alpha.1.1*prt1_1+(1-alpha.0.1)*(1-prt1_1))/(1-prw1_1)
  
  deriv.beta1.g.w1.1=(exp(beta.g)*(1-alpha.1.1)*prt1_1)/(exp(beta.g)*(1-alpha.1.1)*prt1_1 +alpha.0.1*(1-prt1_1))
  
  deriv.beta1.g.w0.1=(exp(beta.g)*alpha.1.1*prt1_1)/(exp(beta.g)*alpha.1.1*prt1_1+(1-alpha.0.1)*(1-prt1_1))
  
  deriv.beta1.g.1=deriv.beta1.g.w1.1*w+(1-w)*deriv.beta1.g.w0.1
  
  exp.g.1=exp.g.w1.1*w+(1-w)*exp.g.w0.1
  
  ################ conditional probability calculation
  term1.2=exp.g.1*exp(.prSum(beta.z,z)) 
  term2.2=matrix(term1.2, ncol=(capm+1), byrow=T)
  term3.2=rep(apply(term2.2, 1, sum), each=(capm+1))
  cond.prob.2=term1.2/term3.2
  
  #
  mid.term.1=0
  
  ##Here the estimates for the U equations are estimates
  uvec.1=matrix(NA,nrow=n,byrow=T)
  dim(uvec.1)
  
  #U equations for gamma
  #Need terms involving the instrument and the stratification variables
  #Will need for loops because we don't know the number of strat variables and instruments
  term11=(
    (1-y)*(1-alpha.0.1-alpha.1.1)*prt1_1*
      (1-prt1_1)*(w-prw1_1)/(prw1_1*(1-prw1_1))
  )
  
  term1.mat=matrix(term11, ncol=(capm+1), byrow=T)
  #dim(term1.mat)
  #summary(term1.mat)
  #instruments
  uvec.1[,1]=apply(term1.mat,1,sum)
  uvecxsterm=term11*xs
  uvecxs=as.matrix(uvecxsterm,nrow=n*(capm+1),byrow=FALSE)
  uvecxs.vec=as.vector(uvecxs)
  uvecxs.mat=matrix(uvecxs.vec,nrow=n,ncol=(capm+1)*len.xs,byrow=TRUE)
  uvecxs.array=apply(array(uvecxs.mat,dim=c(dim(uvecxs.mat)/c(1,len.xs),len.xs)),MARGIN=c(1,3),FUN=sum)
  uvec.1=cbind(uvec.1,uvecxs.array)
  
  #
  #
  #stratification variables
  uvecsvterm=term11*sv
  uvecsv=as.matrix(uvecsvterm,nrow=n*(capm+1),byrow=FALSE)
  uvecsv.vec=as.vector(uvecsv)
  uvecsv.mat=matrix(uvecsv.vec,nrow=n,ncol=(capm+1)*len.sv,byrow=TRUE)
  uvecsv.array=apply(array(uvecsv.mat,dim=c(dim(uvecsv.mat)/c(1,len.sv),len.sv)),MARGIN=c(1,3),FUN=sum)
  uvec.1=cbind(uvec.1,uvecsv.array)

  #
  #
  #U equations for eta
  term2= (1-y)*(w-prw1_1)/(prw1_1*(1-prw1_1))
  
  uveceta1=term2*( alpha.0.1*(1-alpha.0.1)-alpha.0.1*(1-alpha.0.1-alpha.1.1)*prt1_1 )
  uveceta1.mat=matrix(uveceta1,ncol=(capm+1),byrow=T)
  uveceta1.sum=apply(uveceta1.mat,1,sum)
  uvec.1=cbind(uvec.1,uveceta1.sum)
  uveceta2= term2*( -alpha.0.1*alpha.1.1-alpha.1.1*(1-alpha.0.1-alpha.1.1)*prt1_1) 
  uveceta2.mat=matrix(uveceta2,ncol=(capm+1),byrow=T)
  uveceta2.sum=apply(uveceta2.mat,1,sum)
  uvec.1=cbind(uvec.1,uveceta2.sum)

  
  #U equations for beta
  newterm1.1=y-cond.prob.2
  
  mult.1=(newterm1.1)*deriv.beta1.g.1
  mult.mat.1=matrix(mult.1,ncol=(capm+1),byrow=T)
  mult.mat.1.sum=apply(mult.mat.1,1,sum)
  uvec.1=cbind(uvec.1,mult.mat.1.sum)
  #
  uvecbetazterm=newterm1.1*z
  uvecbetaz=as.matrix(uvecbetazterm,nrow=n*(capm+1),byrow=FALSE)
  uvecbetaz.vec=as.vector(uvecbetaz)
  uvecbetaz.mat=matrix(uvecbetaz.vec,nrow=n,ncol=(capm+1)*len.z,byrow=TRUE)
  uvecbetaz.array=apply(array(uvecbetaz.mat,dim=c(dim(uvecbetaz.mat)/c(1,len.z),len.z)),MARGIN=c(1,3),FUN=sum)
  uvec.1=cbind(uvec.1,uvecbetaz.array)


    ##After calculating the uvec, we can now calculate the middle term 
  #for asymptotic variance
  ##mid.term
  
  mid.term1.1=t(uvec.1)%*%uvec.1
  mid.term1.1=mid.term1.1/n
  mid.term.1=mid.term1.1
  
  return(list(uvec.1=uvec.1,mid.term.1=mid.term.1))
}
