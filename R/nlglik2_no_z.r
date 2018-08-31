#for testing purposes
#no.par.two is number of beta parameters, including beta1, 
#and number of beta for each z covariate

#Run this first for testing


#first calculate starting values from naive analysis
.part2noz=function(no.cases, no.controls.per.set, y,sv,w,xs,prt1_1,alpha.0.1,alpha.1.1,prw1_1,naiveestimate,len.xs,len.sv,no.par.one,no.par.two,n,capm,nallpara){
nlglik2=function(beta){
  beta.g=beta[1]
  
  exp.g.w1.1= (exp(beta.g)*(1-alpha.1.1)*prt1_1 +alpha.0.1*(1-prt1_1))/ prw1_1
  exp.g.w0.1= (exp(beta.g)*alpha.1.1*prt1_1+(1-alpha.0.1)*(1-prt1_1))/(1-prw1_1)
  
  deriv.beta1.g.w1.1=(exp(beta.g)*(1-alpha.1.1)*prt1_1)/(exp(beta.g)*(1-alpha.1.1)*prt1_1 +alpha.0.1*(1-prt1_1))
  
  deriv.beta1.g.w0.1=(exp(beta.g)*alpha.1.1*prt1_1)/(exp(beta.g)*alpha.1.1*prt1_1+(1-alpha.0.1)*(1-prt1_1))
  
  deriv.beta1.g.1=deriv.beta1.g.w1.1*w+(1-w)*deriv.beta1.g.w0.1
  
  exp.g.1=exp.g.w1.1*w+(1-w)*exp.g.w0.1
  
  ################ conditional probability calculation
  term1.2=exp.g.1
  term2.2=matrix(term1.2, ncol=(capm+1), byrow=T)
  term3.2=rep(apply(term2.2, 1, sum), each=(capm+1))
  cond.prob.2=term1.2/term3.2
  
  #
  m.2=0*diag(no.par.two)
  newterm1.2=cond.prob.2*deriv.beta1.g.1
  
  
  #flatten the newterm2.2 to create a matrix
  #make a function
  #
  #
  am1.2=matrix(newterm1.2, ncol=(capm+1), byrow=T)
  newterm3.2=apply(am1.2, 1, sum)
  #
  
  #
  #
  #Calculate terms for first row of m*m matrix
  m.2=sum(cond.prob.2*deriv.beta1.g.1^2)-sum(newterm3.2*newterm3.2)
  
  
  ee=-sum(y*log(cond.prob.2))-0.5*log(m.2) 
  return(ee)
  #End function
  }

out3=optim(naiveestimate, nlglik2, method="L-BFGS-B", control=list(maxit=1500))
beta=out3$par

#End of function
return(list(beta=beta,nlglik.val=out3$value,converge=out3$convergence))
}
