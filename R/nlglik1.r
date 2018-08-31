#Begin file
.part1 <- function(no.cases, no.controls.per.set, y,sv,w,xs,startingvalues,len.xs,len.sv,no.par.one,no.par.two,n,capm,nallpara ){
#set variables and variable dimension

##
#For debugging
#
#par=runif(no.par.one, -0.75, 0.75)
###############################################################################
#Begin estimation 
#Specify variables
  s.ad=1+len.xs+len.sv+2
  index1=len.xs+1
  index21=len.xs+2
  index22=s.ad-2
  index23=s.ad-1  

#1st likelihood

# Only 1 xs but more than 1sv
#
#
nlglik_1_1xs=function(par){
  #Extract parameters from par object
  par.int=par[1]
  par.xs=par[2:(1+len.xs)]
  par.sv=par[(2+len.xs):(1+len.xs+len.sv)]
  par.eta=par[(no.par.one-1):(no.par.one)]
  par.set=c(par.int,par.xs,par.sv,par.eta)
  index1=len.xs+1
  index21=len.xs+2
  index22=s.ad-2
  index23=s.ad-1  
  s.ad=1+len.xs+len.sv+2
  #
  ## prpi=pr(X=1|X^*, Z) # a logistic model is assumed 
  prt1_1=1/(1+exp(-par.int+.nrSum(par.xs,xs)+.nrSum(par.sv,sv)))
  #
  ## alpha.0.1=pr(W=1|X=0, Z)     
  alpha.0.1=1/(1+exp(-par.eta[1])+exp(par.eta[2]-par.eta[1]) )
  ##   alpha.1.1=pr(W=0|X=1, Z)
  alpha.1.1=1/(1+exp(-par.eta[2])+exp(par.eta[1]-par.eta[2]) )
  
  
  #
  #In the event that we get an extremely high or extremely low
  #probability of treatment, set it to a specific value
  prt1_1[prt1_1>0.99999]=0.99999
  prt1_1[prt1_1<0.00001]=0.00001
  #
  # prw.1=pr(W=1|X^*, Z)
  prw1_1=alpha.0.1+(1-alpha.0.1-alpha.1.1)*prt1_1
  prw1_1[prw1_1>0.99999]=0.99999
  prw1_1[prw1_1<0.00001]=0.00001
  
  ########################
  term1_1=(
    (1-y)*(1-alpha.0.1-alpha.1.1)^2*prt1_1^2*(1-prt1_1)^2/(prw1_1*(1-prw1_1))
  )
  #
  term2_1=(
    (1-y)*(1-alpha.0.1-alpha.1.1)*prt1_1*(1-prt1_1)/(prw1_1*(1-prw1_1))
  )
  # 
  term3_1=(
    (1-y)/(prw1_1*(1-prw1_1))
  )
  #
  qnty1_1=alpha.0.1*(1-alpha.0.1)+(-alpha.0.1*(1-alpha.0.1)+alpha.0.1*alpha.1.1)*prt1_1
  #
  qnty2_1=-alpha.0.1*alpha.1.1+(alpha.0.1*alpha.1.1-alpha.1.1*(1-alpha.1.1))*prt1_1
  ########################
  ##Calculate firth penalty matrix, ie the analytic derivative
  ##maximum size of analytic derivative matrix is 1+len.xs+len.sv+2
  ##1 corresponds to terms in first row containing solely term1_1
  ##2 corresponds to last two rows, which contain the term3_1
  a.deriv=0*diag(s.ad)
  a.deriv[1, 1]=a.deriv[1, 1]+sum(term1_1)
  #
  #
  ##Calculate elements for first row of analytic derivative matrix
  for (i in 2:(index1)){
    a.deriv[1, i]=a.deriv[1, i]+sum(term1_1*xs[,i-1])
  }
  for (j in (index21):(index22)){
    a.deriv[1, j]=a.deriv[1, j]+sum(term1_1*sv[,j-(len.xs+1)])
  }
  a.deriv[1, (index23)]=a.deriv[1, (index23)]+sum(term2_1*qnty1_1)
  a.deriv[1, (s.ad)]=a.deriv[1, (s.ad)]+sum(term2_1*qnty2_1)
  #
  #
  #
  ##Calculate elements for  rows involving instrument xs
  for (i in 2:(index1)){
    a.deriv[i, i]=a.deriv[i, i]+sum(term1_1*xs[,i-1]*xs[,i-1])
    #Need to also multiply instrument with each stratification variable
    for (j in (index21):(index22)){
      a.deriv[i, j]=a.deriv[i, j]+sum(term1_1*xs[,i-1]*sv[,j-(index1)])
      }
    a.deriv[i, (index23)]=a.deriv[i, (index23)]+sum(term2_1*qnty1_1*xs[,i-1])
    a.deriv[i, s.ad]=a.deriv[i, s.ad]+sum(term2_1*qnty2_1*xs[,i-1])
    }
  #
  #
  #
  ##Calculate elements for  row involving instrument sv
  for (j in (index21):(index22)){
    a.deriv[j, j]=a.deriv[j, j]+sum(term1_1*sv[,j-(index1)]*sv[,j-(index1)])
    ##If number of stratification variables is 1 skip the above line
    a.deriv[j, (index23)]=a.deriv[j, (index23)]+sum(term2_1*qnty1_1*sv[,j-(index1)])
    a.deriv[j, (s.ad)]=a.deriv[j, s.ad]+sum(term2_1*qnty2_1*sv[,j-(index1)])
    #Need to also multiply stratification variable with other stratification variables
    #If number of instruments is 1 skip the above line
    if (j!=(index22)){
      for (l in (j+1):(index22)){
        a.deriv[j, l]=a.deriv[j, l]+sum(term1_1*sv[,j-(index1)]*sv[,(l-(index1))])
      }
    }    
    }
  #   
  #
  #
  ##Calculate elements for last two rows
  a.deriv[(index23), (index23)]=a.deriv[(index23), (index23)]+sum(term3_1*qnty1_1*qnty1_1)
  a.deriv[(index23), s.ad]=a.deriv[(index23), s.ad]+sum(term3_1*qnty1_1*qnty2_1)
  #
  a.deriv[s.ad, s.ad]=a.deriv[s.ad, s.ad]+sum(term3_1*qnty2_1*qnty2_1)
  #
  #
  #
  
  ##Make analytic derivative symmetric matrix
  
  a.deriv[2:s.ad, 1]=a.deriv[1, 2:s.ad]
  a.deriv[s.ad, (index23)]=a.deriv[(index23), s.ad]
  for (k in 2:(index22)){
    a.deriv[(k+1):s.ad, k]=a.deriv[k, (k+1):s.ad]
    }

  #det(a.deriv)
  neglk_1=- sum ((1-y)*(log(prw1_1)*w+(1-w)*log(1-prw1_1)))-0.5*log(det(a.deriv))
  #
  return(neglk_1)}

# Only 1 sv and more than 1 xs
#
#
nlglik_1_1sv=function(par){
  #Extract parameters from par object
  par.int=par[1]
  par.xs=par[2:(1+len.xs)]
  par.sv=par[(2+len.xs):(1+len.xs+len.sv)]
  par.eta=par[(no.par.one-1):(no.par.one)]
  par.set=c(par.int,par.xs,par.sv,par.eta)
  index1=len.xs+1
  index21=len.xs+2
  index22=s.ad-2
  index23=s.ad-1  
  s.ad=1+len.xs+len.sv+2
  #
  ## prpi=pr(X=1|X^*, Z) # a logistic model is assumed 
  prt1_1=1/(1+exp(-par.int+.nrSum(par.xs,xs)+.nrSum(par.sv,sv)))
  #
  ## alpha.0.1=pr(W=1|X=0, Z)     
  alpha.0.1=1/(1+exp(-par.eta[1])+exp(par.eta[2]-par.eta[1]) )
  ##   alpha.1.1=pr(W=0|X=1, Z)
  alpha.1.1=1/(1+exp(-par.eta[2])+exp(par.eta[1]-par.eta[2]) )
  
  
  #
  #In the event that we get an extremely high or extremely low
  #probability of treatment, set it to a specific value
  prt1_1[prt1_1>0.99999]=0.99999
  prt1_1[prt1_1<0.00001]=0.00001
  #
  # prw.1=pr(W=1|X^*, Z)
  prw1_1=alpha.0.1+(1-alpha.0.1-alpha.1.1)*prt1_1
  prw1_1[prw1_1>0.99999]=0.99999
  prw1_1[prw1_1<0.00001]=0.00001
  
  ########################
  term1_1=(
    (1-y)*(1-alpha.0.1-alpha.1.1)^2*prt1_1^2*(1-prt1_1)^2/(prw1_1*(1-prw1_1))
  )
  #
  term2_1=(
    (1-y)*(1-alpha.0.1-alpha.1.1)*prt1_1*(1-prt1_1)/(prw1_1*(1-prw1_1))
  )
  # 
  term3_1=(
    (1-y)/(prw1_1*(1-prw1_1))
  )
  #
  qnty1_1=alpha.0.1*(1-alpha.0.1)+(-alpha.0.1*(1-alpha.0.1)+alpha.0.1*alpha.1.1)*prt1_1
  #
  qnty2_1=-alpha.0.1*alpha.1.1+(alpha.0.1*alpha.1.1-alpha.1.1*(1-alpha.1.1))*prt1_1
  ########################
  ##Calculate firth penalty matrix, ie the analytic derivative
  ##maximum size of analytic derivative matrix is 1+len.xs+len.sv+2
  ##1 corresponds to terms in first row containing solely term1_1
  ##2 corresponds to last two rows, which contain the term3_1
  a.deriv=0*diag(s.ad)
  a.deriv[1, 1]=a.deriv[1, 1]+sum(term1_1)
  #
  #
  ##Calculate elements for first row of analytic derivative matrix
  for (i in 2:(index1)){
    a.deriv[1, i]=a.deriv[1, i]+sum(term1_1*xs[,i-1])
  }
  for (j in (index21):(index22)){
    a.deriv[1, j]=a.deriv[1, j]+sum(term1_1*sv[,j-(len.xs+1)])
  }
  a.deriv[1, (index23)]=a.deriv[1, (index23)]+sum(term2_1*qnty1_1)
  a.deriv[1, (s.ad)]=a.deriv[1, (s.ad)]+sum(term2_1*qnty2_1)
  #
  #
  #
  ##Calculate elements for  rows involving instrument xs
  for (i in 2:(index1)){
    a.deriv[i, i]=a.deriv[i, i]+sum(term1_1*xs[,i-1]*xs[,i-1])
    #Need to also multiply instrument with each stratification variable
    for (j in (index21):(index22)){
      a.deriv[i, j]=a.deriv[i, j]+sum(term1_1*xs[,i-1]*sv[,j-(index1)])
    }
    a.deriv[i, (index23)]=a.deriv[i, (index23)]+sum(term2_1*qnty1_1*xs[,i-1])
    a.deriv[i, s.ad]=a.deriv[i, s.ad]+sum(term2_1*qnty2_1*xs[,i-1])
    #Need to also multiply instrument with other instrument variable
    #If number of instruments is 1 skip the above line
    if (i!=(index1)){
      for (k in (i+1):(index1)){
        a.deriv[i, k]=a.deriv[i, k]+sum(term1_1*xs[,k-1]*xs[,i-1])
      }
    }
  }
  #
  #
  #
  ##Calculate elements for  row involving instrument sv
  for (j in (index21):(index22)){
    a.deriv[j, j]=a.deriv[j, j]+sum(term1_1*sv[,j-(index1)]*sv[,j-(index1)])
    ##If number of stratification variables is 1 skip the above line
    a.deriv[j, (index23)]=a.deriv[j, (index23)]+sum(term2_1*qnty1_1*sv[,j-(index1)])
    a.deriv[j, (s.ad)]=a.deriv[j, s.ad]+sum(term2_1*qnty2_1*sv[,j-(index1)])
  }
  #   
  #
  #
  ##Calculate elements for last two rows
  a.deriv[(index23), (index23)]=a.deriv[(index23), (index23)]+sum(term3_1*qnty1_1*qnty1_1)
  a.deriv[(index23), s.ad]=a.deriv[(index23), s.ad]+sum(term3_1*qnty1_1*qnty2_1)
  #
  a.deriv[s.ad, s.ad]=a.deriv[s.ad, s.ad]+sum(term3_1*qnty2_1*qnty2_1)
  #
  #
  #
  
  ##Make analytic derivative symmetric matrix
  
  a.deriv[2:s.ad, 1]=a.deriv[1, 2:s.ad]
  a.deriv[s.ad, (index23)]=a.deriv[(index23), s.ad]
  for (k in 2:(index22)){
    a.deriv[(k+1):s.ad, k]=a.deriv[k, (k+1):s.ad]
  }
  
  #det(a.deriv)
  neglk_1=- sum ((1-y)*(log(prw1_1)*w+(1-w)*log(1-prw1_1)))-0.5*log(det(a.deriv))
  #
  return(neglk_1)}

# Only 1 xs and 1 sv
#
#
nlglik_1_1xs1sv=function(par){
  #Extract parameters from par object
  par.int=par[1]
  par.xs=par[2:(1+len.xs)]
  par.sv=par[(2+len.xs):(1+len.xs+len.sv)]
  par.eta=par[(no.par.one-1):(no.par.one)]
  par.set=c(par.int,par.xs,par.sv,par.eta)
  index1=len.xs+1
  index21=len.xs+2
  index22=s.ad-2
  index23=s.ad-1  
  s.ad=1+len.xs+len.sv+2
  #
  ## prpi=pr(X=1|X^*, Z) # a logistic model is assumed 
  prt1_1=1/(1+exp(-par.int+.nrSum(par.xs,xs)+.nrSum(par.sv,sv)))
  #
  ## alpha.0.1=pr(W=1|X=0, Z)     
  alpha.0.1=1/(1+exp(-par.eta[1])+exp(par.eta[2]-par.eta[1]) )
  ##   alpha.1.1=pr(W=0|X=1, Z)
  alpha.1.1=1/(1+exp(-par.eta[2])+exp(par.eta[1]-par.eta[2]) )
  
  
  #
  #In the event that we get an extremely high or extremely low
  #probability of treatment, set it to a specific value
  prt1_1[prt1_1>0.99999]=0.99999
  prt1_1[prt1_1<0.00001]=0.00001
  #
  # prw.1=pr(W=1|X^*, Z)
  prw1_1=alpha.0.1+(1-alpha.0.1-alpha.1.1)*prt1_1
  prw1_1[prw1_1>0.99999]=0.99999
  prw1_1[prw1_1<0.00001]=0.00001
  
  ########################
  term1_1=(
    (1-y)*(1-alpha.0.1-alpha.1.1)^2*prt1_1^2*(1-prt1_1)^2/(prw1_1*(1-prw1_1))
  )
  #
  term2_1=(
    (1-y)*(1-alpha.0.1-alpha.1.1)*prt1_1*(1-prt1_1)/(prw1_1*(1-prw1_1))
  )
  # 
  term3_1=(
    (1-y)/(prw1_1*(1-prw1_1))
  )
  #
  qnty1_1=alpha.0.1*(1-alpha.0.1)+(-alpha.0.1*(1-alpha.0.1)+alpha.0.1*alpha.1.1)*prt1_1
  #
  qnty2_1=-alpha.0.1*alpha.1.1+(alpha.0.1*alpha.1.1-alpha.1.1*(1-alpha.1.1))*prt1_1
  ########################
  ##Calculate firth penalty matrix, ie the analytic derivative
  ##maximum size of analytic derivative matrix is 1+len.xs+len.sv+2
  ##1 corresponds to terms in first row containing solely term1_1
  ##2 corresponds to last two rows, which contain the term3_1
  a.deriv=0*diag(s.ad)
  a.deriv[1, 1]=a.deriv[1, 1]+sum(term1_1)
  #
  #
  ##Calculate elements for first row of analytic derivative matrix
  for (i in 2:(index1)){
    a.deriv[1, i]=a.deriv[1, i]+sum(term1_1*xs[,i-1])
  }
  for (j in (index21):(index22)){
    a.deriv[1, j]=a.deriv[1, j]+sum(term1_1*sv[,j-(len.xs+1)])
  }
  a.deriv[1, (index23)]=a.deriv[1, (index23)]+sum(term2_1*qnty1_1)
  a.deriv[1, (s.ad)]=a.deriv[1, (s.ad)]+sum(term2_1*qnty2_1)
  #
  #
  #
  ##Calculate elements for  rows involving instrument xs
  for (i in 2:(index1)){
    a.deriv[i, i]=a.deriv[i, i]+sum(term1_1*xs[,i-1]*xs[,i-1])
    #Need to also multiply instrument with each stratification variable
    for (j in (index21):(index22)){
      a.deriv[i, j]=a.deriv[i, j]+sum(term1_1*xs[,i-1]*sv[,j-(index1)])
    }
    a.deriv[i, (index23)]=a.deriv[i, (index23)]+sum(term2_1*qnty1_1*xs[,i-1])
    a.deriv[i, s.ad]=a.deriv[i, s.ad]+sum(term2_1*qnty2_1*xs[,i-1])
  }
  #
  #
  #
  ##Calculate elements for  row involving instrument sv
  for (j in (index21):(index22)){
    a.deriv[j, j]=a.deriv[j, j]+sum(term1_1*sv[,j-(index1)]*sv[,j-(index1)])
    ##If number of stratification variables is 1 skip the above line
    a.deriv[j, (index23)]=a.deriv[j, (index23)]+sum(term2_1*qnty1_1*sv[,j-(index1)])
    a.deriv[j, (s.ad)]=a.deriv[j, s.ad]+sum(term2_1*qnty2_1*sv[,j-(index1)])
  }
  #   
  #
  #
  ##Calculate elements for last two rows
  a.deriv[(index23), (index23)]=a.deriv[(index23), (index23)]+sum(term3_1*qnty1_1*qnty1_1)
  a.deriv[(index23), s.ad]=a.deriv[(index23), s.ad]+sum(term3_1*qnty1_1*qnty2_1)
  #
  a.deriv[s.ad, s.ad]=a.deriv[s.ad, s.ad]+sum(term3_1*qnty2_1*qnty2_1)
  #
  #
  #
  
  ##Make analytic derivative symmetric matrix
  
  a.deriv[2:s.ad, 1]=a.deriv[1, 2:s.ad]
  a.deriv[s.ad, (index23)]=a.deriv[(index23), s.ad]
  for (k in 2:(index22)){
    a.deriv[(k+1):s.ad, k]=a.deriv[k, (k+1):s.ad]
  }
  
  #det(a.deriv)
  neglk_1=- sum ((1-y)*(log(prw1_1)*w+(1-w)*log(1-prw1_1)))-0.5*log(det(a.deriv))
  #
  return(neglk_1)}

# More than 1 xs and more than 1 sv
#
#
nlglik_1_xssv=function(par){
  #Extract parameters from par object
  par.int=par[1]
  par.xs=par[2:(1+len.xs)]
  par.sv=par[(2+len.xs):(1+len.xs+len.sv)]
  par.eta=par[(no.par.one-1):(no.par.one)]
  par.set=c(par.int,par.xs,par.sv,par.eta)
  index1=len.xs+1
  index21=len.xs+2
  index22=s.ad-2
  index23=s.ad-1  
  s.ad=1+len.xs+len.sv+2
  #
  ## prpi=pr(X=1|X^*, Z) # a logistic model is assumed 
  prt1_1=1/(1+exp(-par.int+.nrSum(par.xs,xs)+.nrSum(par.sv,sv)))
  #
  ## alpha.0.1=pr(W=1|X=0, Z)     
  alpha.0.1=1/(1+exp(-par.eta[1])+exp(par.eta[2]-par.eta[1]) )
  ##   alpha.1.1=pr(W=0|X=1, Z)
  alpha.1.1=1/(1+exp(-par.eta[2])+exp(par.eta[1]-par.eta[2]) )
  
  
  #
  #In the event that we get an extremely high or extremely low
  #probability of treatment, set it to a specific value
  prt1_1[prt1_1>0.99999]=0.99999
  prt1_1[prt1_1<0.00001]=0.00001
  #
  # prw.1=pr(W=1|X^*, Z)
  prw1_1=alpha.0.1+(1-alpha.0.1-alpha.1.1)*prt1_1
  prw1_1[prw1_1>0.99999]=0.99999
  prw1_1[prw1_1<0.00001]=0.00001
  
  ########################
  term1_1=(
    (1-y)*(1-alpha.0.1-alpha.1.1)^2*prt1_1^2*(1-prt1_1)^2/(prw1_1*(1-prw1_1))
  )
  #
  term2_1=(
    (1-y)*(1-alpha.0.1-alpha.1.1)*prt1_1*(1-prt1_1)/(prw1_1*(1-prw1_1))
  )
  # 
  term3_1=(
    (1-y)/(prw1_1*(1-prw1_1))
  )
  #
  qnty1_1=alpha.0.1*(1-alpha.0.1)+(-alpha.0.1*(1-alpha.0.1)+alpha.0.1*alpha.1.1)*prt1_1
  #
  qnty2_1=-alpha.0.1*alpha.1.1+(alpha.0.1*alpha.1.1-alpha.1.1*(1-alpha.1.1))*prt1_1
  ########################
  ##Calculate firth penalty matrix, ie the analytic derivative
  ##maximum size of analytic derivative matrix is 1+len.xs+len.sv+2
  ##1 corresponds to terms in first row containing solely term1_1
  ##2 corresponds to last two rows, which contain the term3_1
  a.deriv=0*diag(s.ad)
  a.deriv[1, 1]=a.deriv[1, 1]+sum(term1_1)
  #
  #
  ##Calculate elements for first row of analytic derivative matrix
  for (i in 2:(index1)){
    a.deriv[1, i]=a.deriv[1, i]+sum(term1_1*xs[,i-1])
  }
  for (j in (index21):(index22)){
    a.deriv[1, j]=a.deriv[1, j]+sum(term1_1*sv[,j-(len.xs+1)])
  }
  a.deriv[1, (index23)]=a.deriv[1, (index23)]+sum(term2_1*qnty1_1)
  a.deriv[1, (s.ad)]=a.deriv[1, (s.ad)]+sum(term2_1*qnty2_1)
  #
  #
  #
  ##Calculate elements for  rows involving instrument xs
  for (i in 2:(index1)){
    a.deriv[i, i]=a.deriv[i, i]+sum(term1_1*xs[,i-1]*xs[,i-1])
    #Need to also multiply instrument with each stratification variable
    for (j in (index21):(index22)){
      a.deriv[i, j]=a.deriv[i, j]+sum(term1_1*xs[,i-1]*sv[,j-(index1)])
    }
    a.deriv[i, (index23)]=a.deriv[i, (index23)]+sum(term2_1*qnty1_1*xs[,i-1])
    a.deriv[i, s.ad]=a.deriv[i, s.ad]+sum(term2_1*qnty2_1*xs[,i-1])
    #Need to also multiply instrument with other instrument variable
    if (i!=(index1)){
      for (k in (i+1):(index1)){
        a.deriv[i, k]=a.deriv[i, k]+sum(term1_1*xs[,k-1]*xs[,i-1])
      }
    }
  }
  #
  #
  #
  ##Calculate elements for  row involving instrument sv
  for (j in (index21):(index22)){
    a.deriv[j, j]=a.deriv[j, j]+sum(term1_1*sv[,j-(index1)]*sv[,j-(index1)])
    ##If number of stratification variables is 1 skip the above line
    a.deriv[j, (index23)]=a.deriv[j, (index23)]+sum(term2_1*qnty1_1*sv[,j-(index1)])
    a.deriv[j, (s.ad)]=a.deriv[j, s.ad]+sum(term2_1*qnty2_1*sv[,j-(index1)])
    #Need to also multiply stratification variable with other stratification variables
    #If number of instruments is 1 skip the above line
    if (j!=(index22)){
      for (l in (j+1):(index22)){
        a.deriv[j, l]=a.deriv[j, l]+sum(term1_1*sv[,j-(index1)]*sv[,(l-(index1))])
      }
    }    
  }
  #   
  #
  #
  ##Calculate elements for last two rows
  a.deriv[(index23), (index23)]=a.deriv[(index23), (index23)]+sum(term3_1*qnty1_1*qnty1_1)
  a.deriv[(index23), s.ad]=a.deriv[(index23), s.ad]+sum(term3_1*qnty1_1*qnty2_1)
  #
  a.deriv[s.ad, s.ad]=a.deriv[s.ad, s.ad]+sum(term3_1*qnty2_1*qnty2_1)
  #
  #
  #
  
  ##Make analytic derivative symmetric matrix
  
  a.deriv[2:s.ad, 1]=a.deriv[1, 2:s.ad]
  a.deriv[s.ad, (index23)]=a.deriv[(index23), s.ad]
  for (k in 2:(index22)){
    a.deriv[(k+1):s.ad, k]=a.deriv[k, (k+1):s.ad]
  }
  
  #det(a.deriv)
  neglk_1=- sum ((1-y)*(log(prw1_1)*w+(1-w)*log(1-prw1_1)))-0.5*log(det(a.deriv))
  #
  return(neglk_1)}


##################################################################
if (len.xs!=1 & len.sv!=1){
  out2=optim(startingvalues, nlglik_1_xssv, method="L-BFGS-B", control=list(maxit=1500), hessian=T)
} else{
  if (len.xs==1 &len.sv!=1){
    out2=optim(startingvalues, nlglik_1_1xs, method="L-BFGS-B", control=list(maxit=1500), hessian=T)
  } else {
    if (len.xs!=1&len.sv==1){
      out2=optim(startingvalues, nlglik_1_1sv, method="L-BFGS-B", control=list(maxit=1500), hessian=T)
    } else{
        out2=optim(startingvalues, nlglik_1_1xs1sv, method="L-BFGS-B", control=list(maxit=1500), hessian=T)
      }
  }
}

#Now store the optimal values as par to calculate alpha.0.1,alpha.1.1,prt1_1,prw1_1

par2=out2$par
par.int=par2[1]
par.xs=par2[2:(1+len.xs)]
par.sv=par2[(2+len.xs):(1+len.xs+len.sv)]
par.eta=par2[(no.par.one-1):(no.par.one)]


prt1_1=1/(1+exp(-par.int-rowSums(par.xs*xs)-rowSums(par.sv*sv)))
## alpha.0.1=pr(W=1|X=0, Z)     
alpha.0.1=1/(1+exp(-par.eta[1])+exp(par.eta[2]-par.eta[1]) )
##   alpha.1.1=pr(W=0|X=1, Z)
alpha.1.1=1/(1+exp(-par.eta[2])+exp(par.eta[1]-par.eta[2]) )

prw1_1=alpha.0.1+(1-alpha.0.1-alpha.1.1)*prt1_1


#End of function
return(list(prt1_1=prt1_1,alpha.0.1=alpha.0.1,alpha.1.1=alpha.1.1,prw1_1=prw1_1,
            par.int=par.int,par.xs=par.xs,par.sv=par.sv,par.eta=par.eta,hess.mat=out2$hessian,nlglik.val=out2$value,converge=out2$convergence))
}
