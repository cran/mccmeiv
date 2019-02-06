#for testing purposes
#no.par.two is number of beta parameters, including beta1, 
#and number of beta for each z covariate


#first calculate starting values from naive analysis
.part2noa0=function(no.cases, no.controls.per.set, y,sv,w,xs,z,prt1_1,alpha.1.1,prw1_1,naiveestimate,len.xs,len.sv,len.z,no.par.one,no.par.two,n,capm,nallpara){
nlglik2=function(beta){
  
  beta.g=beta[1]
  beta.z=beta[2:no.par.two]
  
  exp.g.w1.1= (exp(beta.g)*(1-alpha.1.1)*prt1_1 +0*(1-prt1_1))/ prw1_1
  exp.g.w0.1= (exp(beta.g)*alpha.1.1*prt1_1+(1-0)*(1-prt1_1))/(1-prw1_1)
  
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
  
  cond.prob.2[cond.prob.2>0.99999]=0.99999
  cond.prob.2[cond.prob.2<0.00001]=0.00001
  
  ee=-sum(y*log(cond.prob.2))
  return(ee)
  #End function
  }

out3=optim(naiveestimate, nlglik2, method="L-BFGS-B", control=list(maxit=1500))
beta=out3$par

#End of function
return(list(beta=beta,nlglik.val=out3$value,converge=out3$convergence))
}