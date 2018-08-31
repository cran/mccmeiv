#This function is used to calculate the test statistics for the estimates from the analysis. Returns a 
#table with the parameters, estimates, std errors, 1-alpha/2 

.finalpart=function(par.set,beta,no.par.one,par.sd,var.cov,alpha,param.labels,n.allpar){
  #Set alpha level
  alpha.crit=1-alpha/2
  eta.0.loc=no.par.one-1
  eta.1.loc=no.par.one
  
  #First calculate estimates for alpha_0 and alpha_1 using delta method
  alpha.0=1/(1+exp(-par.set[eta.0.loc])+exp(par.set[eta.1.loc]-par.set[eta.0.loc]) )# alpha.0=pr(W=1|T=0, Z)
  ##   alpha.1=pr(W=0|T=1, Z)
  alpha.1=1/(1+exp(-par.set[eta.1.loc])+exp(par.set[eta.0.loc]-par.set[eta.1.loc]) )# alpha.1=pr(W=0|T=1, Z)
  #Calculate derivatives
  der.alpha.0.b1=-(-exp(-par.set[eta.0.loc])-exp(par.set[eta.1.loc]-par.set[eta.0.loc]))/(1+exp(-par.set[eta.0.loc])+exp(par.set[eta.1.loc]-par.set[eta.0.loc]))^2
  der.alpha.0.b2=-(exp(par.set[eta.1.loc]-par.set[eta.0.loc]))/(1+exp(-par.set[eta.0.loc])+exp(par.set[eta.1.loc]-par.set[eta.0.loc]))^2
  der.alpha.1.b1=-(exp(par.set[eta.0.loc]-par.set[eta.1.loc]))/(1+exp(-par.set[eta.1.loc])+exp(par.set[eta.0.loc]-par.set[eta.1.loc]))^2
  der.alpha.1.b2=-(-exp(-par.set[eta.1.loc])-exp(par.set[eta.0.loc]-par.set[eta.1.loc]))/(1+exp(-par.set[eta.1.loc])+exp(par.set[eta.0.loc]-par.set[eta.1.loc]))^2
  J0=matrix(c(der.alpha.0.b1,der.alpha.0.b2),ncol=2)
  J1=matrix(c(der.alpha.1.b1,der.alpha.1.b2),ncol=2)
  #Calculate V vector
  V=matrix(var.cov[eta.0.loc:eta.1.loc,eta.0.loc:eta.1.loc],ncol=2)
  #Calculate standard error
  var.err.alpha.0=J0%*%V%*%t(J0)
  var.err.alpha.1=J1%*%V%*%t(J1)
  
  alpha.0
  alpha.1
  alpha0.sd=sqrt(var.err.alpha.0)
  alpha1.sd=sqrt(var.err.alpha.1)
  
  #parameters and standard errors
  estimates=round(c(par.set,beta,alpha.0,alpha.1),digits=5)
  se=round(c(par.sd,alpha0.sd,alpha1.sd),digits=5)
  
  #
  #test statistics 
  test.stat.par.set=round(par.set/par.sd[1:(no.par.one)],digits=5)
  test.stat.beta=round(beta/par.sd[(no.par.one+1):(n.allpar)],digits=5)
  test.stat.alpha0=round(alpha.0/alpha0.sd,digits=5)
  test.stat.alpha1=round(alpha.1/alpha1.sd,digits=5)
  
  test.stat=c(test.stat.par.set,test.stat.beta,test.stat.alpha0,test.stat.alpha1)
  
  #
  #Calculate p values
  pval.par.set=round(2*pnorm(-abs(test.stat.par.set)),digits=5)
  pval.beta=round(2*pnorm(-abs(test.stat.beta)),digits=5)
  pval.alpha0=round(2*pnorm(-abs(test.stat.alpha0)),digits=5)
  pval.alpha1=round(2*pnorm(-abs(test.stat.alpha1)),digits=5)
  
  pval.stat=c(pval.par.set,pval.beta,pval.alpha0,pval.alpha1)
  
  #
  #Calculate 1-alpha CI
  crit.val=qnorm(alpha.crit)
  CI.par.set.lower=round(par.set-crit.val*par.sd[1:(no.par.one)],digits=5)
  CI.par.set.upper=round(par.set+crit.val*par.sd[1:(no.par.one)],digits=5)
  
  
  CI.beta.lower=round(beta-crit.val*par.sd[(no.par.one+1):(n.allpar)],digits=5)
  CI.beta.upper=round(beta+crit.val*par.sd[(no.par.one+1):(n.allpar)],digits=5)
  
  CI.alpha0.lower=round(alpha.0-crit.val*alpha0.sd,digits=5)
  CI.alpha0.upper=round(alpha.0+crit.val*alpha0.sd,digits=5)
  
  CI.alpha1.lower=round(alpha.1-crit.val*alpha1.sd,digits=5)
  CI.alpha1.upper=round(alpha.1+crit.val*alpha1.sd,digits=5)
  
  CI.lower=c(CI.par.set.lower,CI.beta.lower,CI.alpha0.lower,CI.alpha1.lower)
  CI.upper=c(CI.par.set.upper,CI.beta.upper,CI.alpha0.upper,CI.alpha1.upper)
  
  #
  #create labels for output
  finalresults.output=cbind(param.labels,estimates, se,as.numeric(test.stat),as.numeric(CI.lower),as.numeric(CI.upper),as.numeric(pval.stat))
  colnames(finalresults.output)=c("Parameter Estimates","Value","S.E.", "Z_Value","CI_Lower_Limit","CI_Upper_Limit","P-Value")

  return(result=finalresults.output)

}
