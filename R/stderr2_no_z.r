####### Derivative calculation using numerical approach #####
##After calculating the middle term, we now need to calculate A hat
##Once we have A hat we then  now calculate asymptotic variance
################
.part4noz=function(mid.term.1,hess.mat,y,sv,w,xs,allpara,len.xs,len.sv,no.par.one,no.par.two,n,capm,nallpara){
  
  forderiveq2=function(allpara){
  #Convert data into vectors
   
  #extract parameters
  par=allpara[1:no.par.one]
  beta= allpara[(no.par.one+1):nallpara]
  
  #extract parameters
  par.int=par[1]
  par.xs=par[2:(1+len.xs)]
  par.sv=par[(2+len.xs):(1+len.xs+len.sv)]
  par.eta=par[(no.par.one-1):(no.par.one)]

  #
  ## prpi=pr(X=1|X^*, Z) # a logistic model is assumed 
  prt1_1=1/(1+exp(-par.int+.nrSum(par.xs,xs)+.nrSum(par.sv,sv)))
  #
  ## alpha.0.1=pr(W=1|X=0, Z)     
  alpha.0.1=1/(1+exp(-par.eta[1])+exp(par.eta[2]-par.eta[1]) )
  ##   alpha.1.1=pr(W=0|X=1, Z)
  alpha.1.1=1/(1+exp(-par.eta[2])+exp(par.eta[1]-par.eta[2]) )
  
  #
  # prw.1=pr(W=1|X^*, Z)
  prw1_1=alpha.0.1+(1-alpha.0.1-alpha.1.1)*prt1_1
  
  #Function used to find optimal values of beta


    beta.g=beta[1]


    #calculate terms for derivatives
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
  
#
#
#
deriv.1=0*diag(nallpara)#jacobian(forderiv, allpara)/n
deriv.1[1:no.par.one, 1:no.par.one]=-hess.mat/n
out5.1=hessian(forderiveq2, allpara)
deriv.1[(no.par.one+1):nallpara, ]=-out5.1[(no.par.one+1):nallpara, ]/n

##########
##########
##Finally calculate the A hat, then calculate the asymptotic variance
inv.deriv.1=solve(deriv.1)
var.cov.1= inv.deriv.1%*%mid.term.1%*%t(inv.deriv.1)/n
theta.sd.1=sqrt(diag(var.cov.1))

return(list(theta.sd.1=theta.sd.1,var.cov.1=var.cov.1,deriv.1=deriv.1))

#End function
  }
