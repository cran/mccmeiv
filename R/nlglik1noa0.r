#Begin file
.part1noa0 <- function(no.cases, no.controls.per.set, y,sv,w,xs,z,startingvalues,len.xs,len.sv,len.z,no.par.one,no.par.two,n,capm,nallpara ){
#set variables and variable dimension

##
#For debugging
#
#par=runif(no.par.one, -0.75, 0.75)
###############################################################################
#Begin estimation 

#1st likelihood
#
nlglik_1=function(par){
  #Extract parameters from par object
  par.int=par[1]
  par.xs=par[2:(1+len.xs)]
  par.sv=par[(2+len.xs):(1+len.xs+len.sv)]
  par.z=par[(2+len.xs+len.sv):(1+len.xs+len.sv+len.z)]
  par.eta=par[no.par.one]
  par.set=c(par.int,par.xs,par.sv,par.z,par.eta)
  #
  ## prpi=pr(X=1|X^*, Z) # a logistic model is assumed 
  #prt1_1=1/(1+exp(-par.int+nrSum(par.xs,xs)+nrSum(par.sv,sv)+nrSum(par.z,z)))
  sum.xs=as.vector(as.matrix(xs)%*%par.xs)
  sum.sv=as.vector(as.matrix(sv)%*%par.sv)
  sum.z=as.vector(as.matrix(z)%*%par.z)
  prt1_1=1/(1+exp(-par.int-sum.xs-sum.sv-sum.z))
  #
  ## alpha.0.1=pr(W=1|X=0, Z)=0
  ##   alpha.1.1=pr(W=0|X=1, Z)
  alpha.1.1=1/(1+exp(-par.eta) )
  
  
  #
  #In the event that we get an extremely high or extremely low
  #probability of treatment, set it to a specific value
  prt1_1[prt1_1>0.99999]=0.99999
  prt1_1[prt1_1<0.00001]=0.00001
  #
  # prw.1=pr(W=1|X^*, Z)
  prw1_1=(1-alpha.1.1)*prt1_1
  prw1_1[prw1_1>0.99999]=0.99999
  prw1_1[prw1_1<0.00001]=0.00001
  
  ########################
 

  #det(a.deriv)
  neglk_1=- sum ((1-y)*(log(prw1_1)*w+(1-w)*log(1-prw1_1)))
  #
  return(neglk_1)}

out2=optim(startingvalues, nlglik_1,method="L-BFGS-B", control=list(maxit=1500))
#Now store the optimal values as par to calculate alpha.0.1,alpha.1.1,prt1_1,prw1_1

par2=out2$par
par.int=par2[1]
par.xs=par2[2:(1+len.xs)]
par.sv=par2[(2+len.xs):(1+len.xs+len.sv)]
par.z=par2[(2+len.xs+len.sv):(1+len.xs+len.sv+len.z)]
par.eta=par2[no.par.one]


sum.xs=as.vector(as.matrix(xs)%*%par.xs)
sum.sv=as.vector(as.matrix(sv)%*%par.sv)
sum.z=as.vector(as.matrix(z)%*%par.z)
prt1_1=1/(1+exp(-par.int-sum.xs-sum.sv-sum.z))
## alpha.0.1=pr(W=1|X=0, Z)=0    
##   alpha.1.1=pr(W=0|X=1, Z)
alpha.1.1=1/(1+exp(-par.eta) )

prw1_1=(1-alpha.1.1)*prt1_1


#End of function
return(list(prt1_1=prt1_1,alpha.1.1=alpha.1.1,prw1_1=prw1_1,pars2=par2,
            par.int=par.int,par.xs=par.xs,par.sv=par.sv,par.z=par.z,par.eta=par.eta,
            hess.mat=out2$hessian,nlglik.val=out2$value,converge=out2$convergence))
}
