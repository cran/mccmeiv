#This function is used to calculate the test statistics for the estimates from the analysis. Returns a 
#table with the parameters, estimates, std errors, 1-alpha/2 

.finalpartnoa1=function(par.set,beta,no.par.one,par.sd,var.cov,alpha,param.labels,n.allpar){
  #Set alpha level
  alpha.crit=1-alpha/2
  eta.0.loc=no.par.one
  
  
  #First calculate estimates for alpha_0 and alpha_1 using delta method
  alpha.0=1/(1+exp(-par.set[eta.0.loc]))# alpha.0=pr(W=1|T=0, Z)
  ##   alpha.1=pr(W=0|T=1, Z)=0
  #Calculate derivatives
  der.alpha.0.b1=alpha.0*(1-alpha.0)
  #Calculate V vector
  V=var.cov[eta.0.loc,eta.0.loc]
  #Calculate standard error
  var.err.alpha.0=V*der.alpha.0.b1

  alpha.0
  alpha0.sd=sqrt(var.err.alpha.0)

  #parameters and standard errors
  #estimates=round(c(par.set,beta,alpha.0),digits=5)
  estimates=round(beta,digits=5)
  exp_estimates=round(exp(beta),digits=5)
  #se=round(c(par.sd,alpha0.sd),digits=5)
  se=round(par.sd[(no.par.one+1):(n.allpar)],digits=5)
  #
  #test statistics 
  test.stat.par.set=round(par.set/par.sd[1:(no.par.one)],digits=5)
  test.stat.beta=round(beta/par.sd[(no.par.one+1):(n.allpar)],digits=5)
  test.stat.alpha0=round(alpha.0/alpha0.sd,digits=5)

  #test.stat=c(test.stat.par.set,test.stat.beta,test.stat.alpha0)
  test.stat=c(test.stat.beta)
  
  #
  #Calculate p values
  pval.par.set=round(2*pnorm(-abs(test.stat.par.set)),digits=5)
  pval.beta=round(2*pnorm(-abs(test.stat.beta)),digits=5)
  pval.alpha0=round(2*pnorm(-abs(test.stat.alpha0)),digits=5)

  #pval.stat=c(pval.par.set,pval.beta,pval.alpha0)
  pval.stat=c(pval.beta)
  
  #
  #Calculate 1-alpha CI
  crit.val=qnorm(alpha.crit)
  CI.par.set.lower=round(par.set-crit.val*par.sd[1:(no.par.one)],digits=5)
  CI.par.set.upper=round(par.set+crit.val*par.sd[1:(no.par.one)],digits=5)
  
  
  CI.beta.lower=round(exp(beta-crit.val*par.sd[(no.par.one+1):(n.allpar)]),digits=5)
  CI.beta.upper=round(exp(beta+crit.val*par.sd[(no.par.one+1):(n.allpar)]),digits=5)
  
  CI.alpha0.lower=round(alpha.0-crit.val*alpha0.sd,digits=5)
  CI.alpha0.upper=round(alpha.0+crit.val*alpha0.sd,digits=5)

  #CI.lower=c(CI.par.set.lower,CI.beta.lower,CI.alpha0.lower)
  #CI.upper=c(CI.par.set.upper,CI.beta.upper,CI.alpha0.upper)
  
  CI.lower=c(CI.beta.lower)
  CI.upper=c(CI.beta.upper)
  #
  #create labels for output
  finalresults.output=cbind(as.numeric(estimates), as.numeric(exp_estimates), as.numeric(se),as.numeric(test.stat),as.numeric(pval.stat),as.numeric(CI.lower),as.numeric(CI.upper))
  colnames(finalresults.output)=c("coef","exp(coef)","se(coef)", "z","Pr(>|z|)","lower_ci","upper_ci")
  rownames(finalresults.output)=param.labels
  return(result=finalresults.output)

}