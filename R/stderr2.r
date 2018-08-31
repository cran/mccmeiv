####### Derivative calculation using numerical approach #####
##After calculating the middle term, we now need to calculate A hat
##Once we have A hat we then  now calculate asymptotic variance
################
.part4=function(mid.term.1,hess.mat,y,sv,w,xs,z,allpara,len.xs,len.sv,len.z,no.par.one,no.par.two,n,capm,nallpara){
  z.check=as.numeric(len.z>1)
  forderiveq21=function(allpara){
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
    beta.z=beta[2:no.par.two]

    #calculate terms for derivatives
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
    m.2=0*diag(no.par.two)
    newterm1.2=cond.prob.2*deriv.beta1.g.1
    newterm2.2=cond.prob.2*z
    
    #flatten the newterm2.2 to create a matrix
    #make a function
    
    newterm2.2.flat=as.vector(t(newterm2.2))
    
    #
    #
    am1.2=matrix(newterm1.2, ncol=(capm+1), byrow=T)
    am2.2=matrix(newterm2.2.flat, ncol=((capm+1)*(no.par.two-1)), byrow=T)
    newterm3.2=apply(am1.2, 1, sum)
    #
    #newterm4.2.mat represents the crossterms for the number of covariates z
    newterm4.2.mat=matrix(,nrow=(length(w)/(capm+1)),ncol=(no.par.two-1))
    newterm4.2.mat=data.frame(newterm4.2.mat)
    #When there is more than 1 covariate z, then we want to sum only the columns 
    #associated with each z covariate. This means we would have to skip every 1 or two columns.
    #hence why we use the special sum function.
    #if there is only 1 z covariate, then just use regular apply 
    newterm4.2.mat[,1] =apply(am2.2, 1, sum)
    
    #
    #
    #Calculate terms for first row of m*m matrix
    m.2[1,1]=sum(cond.prob.2*deriv.beta1.g.1^2)-sum(newterm3.2*newterm3.2)
    for (i in 2:no.par.two){
      m.2[1, i]=m.2[1, i]+sum(cond.prob.2*deriv.beta1.g.1*z[,(i-1)])- 
        sum(newterm3.2*newterm4.2.mat[,(i-1)])
    }
    #
    #Calculate terms for rows involving covariates z
    #Calculate diagonal of information matrix 
    m.2[2, 2]=sum(cond.prob.2*z[,1]^2)-sum(newterm4.2.mat[,1]*newterm4.2.mat[,1])
    
    ##Make the matrix symmetric
    m.2[2:no.par.two, 1]=m.2[1, 2:no.par.two]
    m.2[no.par.two, (no.par.two-1)]=m.2[(no.par.two-1), no.par.two]
    
    ee=-sum(y*log(cond.prob.2)) -0.5*log(det(m.2)) 

    return(ee)
    #End function
  }
  
  forderiveq22=function(allpara){
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
    beta.z=beta[2:no.par.two]
    
    #calculate terms for derivatives
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
    m.2=0*diag(no.par.two)
    newterm1.2=cond.prob.2*deriv.beta1.g.1
    newterm2.2=cond.prob.2*z
    
    #flatten the newterm2.2 to create a matrix
    #make a function
    
    newterm2.2.flat=as.vector(t(newterm2.2))
    
    #
    #
    am1.2=matrix(newterm1.2, ncol=(capm+1), byrow=T)
    am2.2=matrix(newterm2.2.flat, ncol=((capm+1)*(no.par.two-1)), byrow=T)
    newterm3.2=apply(am1.2, 1, sum)
    #
    #newterm4.2.mat represents the crossterms for the number of covariates z
    newterm4.2.mat=matrix(,nrow=(length(w)/(capm+1)),ncol=(no.par.two-1))
    newterm4.2.mat=data.frame(newterm4.2.mat)
    #When there is more than 1 covariate z, then we want to sum only the columns 
    #associated with each z covariate. This means we would have to skip every 1 or two columns.
    #hence why we use the special sum function.
    for (i in 1:(no.par.two-1)){
      rowtotal=apply(am2.2,1,function(x) sum(x[seq(i,length(x),(len.z))]))
      newterm4.2.mat[,i]=rowtotal
    } 
    
    #
    #
    #Calculate terms for first row of m*m matrix
    m.2[1,1]=sum(cond.prob.2*deriv.beta1.g.1^2)-sum(newterm3.2*newterm3.2)
    for (i in 2:no.par.two){
      m.2[1, i]=m.2[1, i]+sum(cond.prob.2*deriv.beta1.g.1*z[,(i-1)])- 
        sum(newterm3.2*newterm4.2.mat[,(i-1)])
    }
    #
    #Calculate terms for rows involving covariates z
    for (i in 2:no.par.two){
      #Calculate diagonal of information matrix 
      m.2[i, i]=sum(cond.prob.2*z[,(i-1)]^2)-sum(newterm4.2.mat[,(i-1)]*newterm4.2.mat[,(i-1)])
      #Calculate cross product terms 
      if (i!=no.par.two){
        for (k in (i+1):(no.par.two)){
          m.2[i, k]=sum(cond.prob.2*z[,(k-2)]*z[,(k-1)])-sum(newterm4.2.mat[,(k-2)]*newterm4.2.mat[,(k-1)])
        }  
      }
    }
    ##Make the matrix symmetric
    m.2[2:no.par.two, 1]=m.2[1, 2:no.par.two]
    m.2[no.par.two, (no.par.two-1)]=m.2[(no.par.two-1), no.par.two]
    for (k in 2:(no.par.two-2)){
      m.2[(k+1):no.par.two, k]=m.2[k, (k+1):no.par.two]
    }
    
    ee=-sum(y*log(cond.prob.2)) -0.5*log(det(m.2)) 
    return(ee)
    #End function
  }
  
#
#
#
deriv.1=0*diag(nallpara)#jacobian(forderiv, allpara)/n
deriv.1[1:no.par.one, 1:no.par.one]=-hess.mat/n
if (z.check==1){
  out5.1=hessian(forderiveq22, allpara)
} else {
  out5.1=hessian(forderiveq21, allpara)
  }


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
