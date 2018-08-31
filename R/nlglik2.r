#for testing purposes
#no.par.two is number of beta parameters, including beta1, 
#and number of beta for each z covariate


#first calculate starting values from naive analysis
.part2=function(no.cases, no.controls.per.set, y,sv,w,xs,z,prt1_1,alpha.0.1,alpha.1.1,prw1_1,naiveestimate,len.xs,len.sv,len.z,no.par.one,no.par.two,n,capm,nallpara){
  
  #Select nlglik function depending on whether there is 1 z or more
  z.check=as.numeric(len.z>1)
  
  #
  #nlglik21 : only 1 z covariate
  nlglik21=function(beta){
  
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

  #More than 1 z covariate
  nlglik22=function(beta){
    
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
  
  #Get estimates of beta
  if (z.check==1) {
    out3=optim(naiveestimate, nlglik22, method="L-BFGS-B", control=list(maxit=1500))
  } else {
      out3=optim(naiveestimate, nlglik21, method="L-BFGS-B", control=list(maxit=1500))
    }  
beta=out3$par

#End of function
return(list(beta=beta,nlglik.val=out3$value,converge=out3$convergence))
}
