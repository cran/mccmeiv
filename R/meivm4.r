meivm4 <-
function(y, sv, xs, w, z,sv.factor=NULL,xs.factor=NULL,z.factor=NULL, alpha=0.05,scale=TRUE,setalpha0.to.0=FALSE,setalpha1.to.0=FALSE){
  
  #
  #Check 1: Verify variables have been specified
  if(missing(y)) stop("Stopping program...need to specify response y")  
  if(missing(sv)) stop("Stopping program...need to specify confounder variables sv")  
  if(missing(xs)) stop("Stopping program...need to specify instrument variables xs")
  if(missing(z)) stop("Stopping program...need to specify prognostic factor z")  
  if(missing(w)) stop("Stopping program...need to specify mismeasured variable w")  
  
  #Check 2: Ensure that alpha 0 and alpha 1 are not both set to 0
  if(setalpha0.to.0==TRUE & setalpha1.to.0==TRUE) stop("Stopping program...can only set either alpha 1 or alpha 0 to 0")    

  #Third read in data 
  n=sum(y)
  if(!is.null(dim(w))){
    w=w[,1]
  }
  row.w=length(w)
  #cat("The number of mismeasured observations  is",row.w,"\n")
  row.y=length(y) #Assume y is a vector
  #cat("The number of response observations  is",row.y,"\n")
  
  ####Confounding/stratification variables information
  if (!is.null(sv.factor)){
    #cat("Stratification variables", sv.factor,"are specified as categorical","\n")
    #First scale non factor confounders
    if (scale==TRUE){
      sv.names=colnames(sv)
      if (length(sv.factor)!=length(intersect(sv.factor,sv.names))) stop("Stopping program, a specified stratification factor
                                                                         variable not found.")
      #Determine which variables are not factors
      sv.transform.var=setdiff(sv.names,sv.factor)
      
      #Determine which of the variables that are not declared factors are binary
      binary.names=names(which(lapply(sv[sv.transform.var], function(x) c(length(unique(x)))==2)==TRUE))
      
      final.sv.transform.var=setdiff(sv.transform.var,binary.names)
      
      #cat("Transforming the following stratification variables",sv.transform.var,"\n")
      sv[final.sv.transform.var]= sapply(sv[final.sv.transform.var], function(x) c(scale(x)))
      
      #Include the names of binary variables that will be included in dummy variable creation
      sv.factor=c(sv.factor,binary.names)
    }
    sv[sv.factor]=lapply(sv[sv.factor], factor)
    #Determine which variables are factors
    list.of.sv.factors=names(which(sapply(sv,class)=="factor"))
    #
    #Create new factor variables
    for (i in 1:length(list.of.sv.factors)){
      #
      #Find levels of variable
      #For now assuming reference level is smallest number
      old.levels=sort(unique(as.numeric(do.call(paste0, sv[list.of.sv.factors[i]]))))
      #cat("Old levels of factor variable",list.of.sv.factors[i] ,"are",old.levels,"\n")
      
      #
      #Make reference category of zero by subtracting smallest level of factor 
      #variable from each observation
      new.levels=as.numeric(old.levels)-min(old.levels)
      cat("New levels of factor variable are",new.levels,"\n")
      new.levels.values=as.numeric(do.call(paste0, sv[list.of.sv.factors[i]]))-min(old.levels)
      new.dum.vars=new.levels[2:length(new.levels)]
      #Create new dummy variable
      for (j in 1:length(new.dum.vars)){
        #
        #Create name of new dummy variable
        newdummy=paste(list.of.sv.factors[i],new.dum.vars[j],sep=".")
        #
        #Assign value of new variable 
        sv=cbind(sv,assign(newdummy,as.numeric(new.levels.values==new.levels[(j+1)]))) #create new variable based on dummy
        names(sv)[dim(sv)[2]]=newdummy
      }
      
      #
      #Drop the variable from the dataframe
      sv[list.of.sv.factors[i]]=NULL
      
    }
  } else{
    #cat("No specified stratification factor variables.")
    if (scale==TRUE){
      sv.names=colnames(sv)
      
      binary.names=names(which(lapply(sv[sv.names], function(x) c(length(unique(x)))==2)==TRUE))
      if (!is.null(binary.names)){
        
        final.sv.transform.var=setdiff(sv.names,binary.names)
        
        #Create old labels
        
        #cat("Transforming the following instrumental variables",xs.transform.var,"\n")
        sv[final.sv.transform.var]= lapply(sv[final.sv.transform.var], function(x) c(scale(x)))
        
        #Create dummy variables for binary variables
        sv.factor=binary.names
        
        sv[sv.factor]=lapply(sv[sv.factor], factor)
        
        #Determine which variables are factors
        list.of.sv.factors=names(which(sapply(sv,class)=="factor"))
        #
        #Create new factor variables
        for (i in 1:length(list.of.sv.factors)){
          #
          #Find levels of variable
          #For now assuming reference level is smallest number
          old.levels=sort(unique(as.numeric(do.call(paste0, sv[list.of.sv.factors[i]]))))
          #cat("Old levels of factor variable",list.of.sv.factors[i] ,"are",old.levels,"\n")
          
          #
          #Make reference category of zero by subtracting smallest level of factor 
          #variable from each observation
          new.levels=as.numeric(old.levels)-min(old.levels)
          cat("New levels of factor variable are",new.levels,"\n")
          new.levels.values=as.numeric(do.call(paste0, sv[list.of.sv.factors[i]]))-min(old.levels)
          new.dum.vars=new.levels[2:length(new.levels)]
          #Create new dummy variable
          for (j in 1:length(new.dum.vars)){
            #
            #Create name of new dummy variable
            newdummy=paste(list.of.sv.factors[i],new.dum.vars[j],sep=".")
            #
            #Assign value of new variable 
            sv=cbind(sv,assign(newdummy,as.numeric(new.levels.values==new.levels[(j+1)]))) #create new variable based on dummy
            names(sv)[dim(sv)[2]]=newdummy
          }
          
          #
          #Drop the variable from the dataframe
          sv[list.of.sv.factors[i]]=NULL
          
        }}
      
      #cat("Transforming all stratification variables","\n")
      else{sv=scale(sv)
      }}}
  sv=data.frame(sv) 
  dimsv=dim(sv)
  len.sv=dim(sv)[2]
  if (names(sv)[1]=="X1") {
    names.sv=rep(NA,len.sv)
    for (i in 1:len.sv){
      names.sv[i]=paste("sv",i,sep="")
      colnames(sv)=names.sv
    }
  } else {
    names.sv=names(sv) 
  }
  row.sv=dim(sv)[1]
  #cat("The number of observations in the confounding variables are",row.sv,"\n")
  
  #
  #
  #Instrument variables information
  if (!is.null(xs.factor)){
    #cat("Instrumental variables", xs.factor,"are specified as categorical")
    #First scale non factor confounders
    if (scale==TRUE){
      xs.names=colnames(xs)
      if (length(xs.factor)!=length(intersect(xs.factor,xs.names))) stop("Stopping program, a specified instrumental factor 
                                                                         variable not found.")
      xs.transform.var=setdiff(xs.names,xs.factor)
      
      #Determine which of the variables that are not declared factors are binary
      binary.names=names(which(lapply(xs[xs.transform.var], function(x) c(length(unique(x)))==2)==TRUE))
      
      final.xs.transform.var=setdiff(xs.transform.var,binary.names)
      #cat("Transforming the following instrumental variables",xs.transform.var,"\n")
      xs[final.xs.transform.var]= lapply(xs[final.xs.transform.var], function(x) c(scale(x)))
      
      #Include the names of binary variables that will be included in dummy variable creation
      xs.factor=c(xs.factor,binary.names)
    }
    xs[xs.factor]=lapply(xs[xs.factor], factor)
    #Determine which variables are factors
    list.of.xs.factors=names(which(sapply(xs,class)=="factor"))
    #
    #Create new factor variables
    for (i in 1:length(list.of.xs.factors)){
      #
      #Find levels of variable
      #For now assuming reference level is smallest number
      old.levels=sort(unique(as.numeric(do.call(paste0, xs[list.of.xs.factors[i]]))))
      #cat("Old levels of factor variable",list.of.xs.factors[i] ,"are",old.levels,"\n")
      
      #
      #Make reference category of zero by subtracting smallest level of factor 
      #variable from each observation
      new.levels=as.numeric(old.levels)-min(old.levels)
      #cat("New levels of factor variable are",new.levels,"\n")
      new.levels.values=as.numeric(do.call(paste0, xs[list.of.xs.factors[i]]))-min(old.levels)
      new.dum.vars=new.levels[2:length(new.levels)]
      #Create new dummy variable
      for (j in 1:length(new.dum.vars)){
        #
        #Create name of new dummy variable
        newdummy=paste(list.of.xs.factors[i],new.dum.vars[j],sep=".")
        #
        #Assign value of new variable 
        xs=cbind(xs,assign(newdummy,as.numeric(new.levels.values==new.levels[(j+1)]))) #create new variable based on dummy
        names(xs)[dim(xs)[2]]=newdummy
      }
      
      #
      #Drop the variable from the dataframe
      xs[list.of.xs.factors[i]]=NULL
      
    }
  } else{
    #cat("No specified instrumental factor variables.")
    if (scale==TRUE){
      xs.names=colnames(xs)
      
      binary.names=names(which(lapply(xs[xs.names], function(x) c(length(unique(x)))==2)==TRUE))
      if (!is.null(binary.names)){
        
        final.xs.transform.var=setdiff(xs.names,binary.names)
        
        #Create old labels
        
        #cat("Transforming the following instrumental variables",xs.transform.var,"\n")
        xs[final.xs.transform.var]= lapply(xs[final.xs.transform.var], function(x) c(scale(x)))
        
        #Create dummy variables for binary variables
        xs.factor=binary.names
        
        xs[xs.factor]=lapply(xs[xs.factor], factor)
        #Determine which variables are factors
        list.of.xs.factors=names(which(sapply(xs,class)=="factor"))
        #
        #Create new factor variables
        for (i in 1:length(list.of.xs.factors)){
          #
          #Find levels of variable
          #For now assuming reference level is smallest number
          old.levels=sort(unique(as.numeric(do.call(paste0, xs[list.of.xs.factors[i]]))))
          #cat("Old levels of factor variable",list.of.xs.factors[i] ,"are",old.levels,"\n")
          
          #
          #Make reference category of zero by subtracting smallest level of factor 
          #variable from each observation
          new.levels=as.numeric(old.levels)-min(old.levels)
          #cat("New levels of factor variable are",new.levels,"\n")
          new.levels.values=as.numeric(do.call(paste0, xs[list.of.xs.factors[i]]))-min(old.levels)
          new.dum.vars=new.levels[2:length(new.levels)]
          #Create new dummy variable
          for (j in 1:length(new.dum.vars)){
            #
            #Create name of new dummy variable
            newdummy=paste(list.of.xs.factors[i],new.dum.vars[j],sep=".")
            #
            #Assign value of new variable 
            xs=cbind(xs,assign(newdummy,as.numeric(new.levels.values==new.levels[(j+1)]))) #create new variable based on dummy
            names(xs)[dim(xs)[2]]=newdummy
          }
          
          #
          #Drop the variable from the dataframe
          xs[list.of.xs.factors[i]]=NULL
        }}
      else{
        
        #cat("Transforming all instrumental variables","\n")
        xs=scale(xs)}
    }}
  xs=data.frame(xs) 
  dimxs=dim(xs)
  len.xs=dim(xs)[2]
  row.xs=dim(xs)[1]
  #cat("The number of observations in instruments variables are",row.xs,"\n")
  if (names(xs)[1]=="X1") {
    names.xs=rep(NA,len.xs)
    for (i in 1:len.xs){
      names.xs[i]=paste("xs",i,sep="")
      colnames(xs)=names.xs
    }
  } else {
    names.xs=names(xs)
  }
  
  #
  #
  #Miscellaneous covariate variables information
  if (!is.null(z.factor)){
    #cat("Prognostic factors z", z.factor,"are specified as categorical","\n")
    #First scale non factor confounders
    if (scale==TRUE){
      z.names=colnames(z)
      if (length(z.factor)!=length(intersect(z.factor,z.names))) stop("Stopping program, a specified prognostic factor 
                                                                      variable not found.")
      z.transform.var=setdiff(z.names,z.factor)
      
      #Determine which of the variables that are not declared factors are binary
      binary.names=names(which(lapply(z[z.transform.var], function(x) c(length(unique(x)))==2)==TRUE))
      
      final.z.transform.var=setdiff(z.transform.var,binary.names)
      
      z[final.z.transform.var]= lapply(z[final.z.transform.var], function(x) c(scale(x)))
      
      z.factor=c(z.factor,binary.names)
    }
    z[z.factor]=lapply(z[z.factor], factor)
    #Determine which variables are factors
    list.of.z.factors=names(which(sapply(z,class)=="factor"))
    #
    #Create new factor variables
    for (i in 1:length(list.of.z.factors)){
      #
      #Find levels of variable
      #For now assuming reference level is smallest number
      old.levels=sort(unique(as.numeric(do.call(paste0, z[list.of.z.factors[i]]))))
      #cat("Old levels of factor variable",list.of.z.factors[i] ,"are",old.levels,"\n")
      
      #
      #Make reference category of zero by subtracting smallest level of factor 
      #variable from each observation
      new.levels=as.numeric(old.levels)-min(old.levels)
      #cat("New levels of factor variable are",new.levels,"\n")
      new.levels.values=as.numeric(do.call(paste0, z[list.of.z.factors[i]]))-min(old.levels)
      new.dum.vars=new.levels[2:length(new.levels)]
      #Create new dummy variable
      for (j in 1:length(new.dum.vars)){
        #
        #Create name of new dummy variable
        newdummy=paste(list.of.z.factors[i],new.dum.vars[j],sep=".")
        #
        #Assign value of new variable 
        z=cbind(z,assign(newdummy,as.numeric(new.levels.values==new.levels[(j+1)]))) #create new variable based on dummy
        names(z)[dim(z)[2]]=newdummy
      }
      
      #
      #Drop the variable from the dataframe
      z[list.of.z.factors[i]]=NULL
      
    }
  } else{
    #cat("No specified z covariate factor variables.","\n")
    if (scale==TRUE){
      z.names=colnames(z)
      
      binary.names=names(which(lapply(z[z.names], function(x) c(length(unique(x)))==2)==TRUE))
      
      if (!is.null(binary.names)){
        
        final.z.transform.var=setdiff(z.names,binary.names)
        
        #Create old labels
        
        #cat("Transforming the following instrumental variables",xs.transform.var,"\n")
        z[final.z.transform.var]= lapply(z[final.z.transform.var], function(x) c(scale(x)))
        
        #Create dummy variables for binary variables
        z.factor=binary.names
        
        z[z.factor]=lapply(z[z.factor], factor)
        #Determine which variables are factors
        list.of.z.factors=names(which(sapply(z,class)=="factor"))
        #
        #Create new factor variables
        for (i in 1:length(list.of.xs.factors)){
          #
          #Find levels of variable
          #For now assuming reference level is smallest number
          old.levels=sort(unique(as.numeric(do.call(paste0, xs[list.of.xs.factors[i]]))))
          #cat("Old levels of factor variable",list.of.xs.factors[i] ,"are",old.levels,"\n")
          
          #
          #Make reference category of zero by subtracting smallest level of factor 
          #variable from each observation
          new.levels=as.numeric(old.levels)-min(old.levels)
          #cat("New levels of factor variable are",new.levels,"\n")
          new.levels.values=as.numeric(do.call(paste0, z[list.of.z.factors[i]]))-min(old.levels)
          new.dum.vars=new.levels[2:length(new.levels)]
          #Create new dummy variable
          for (j in 1:length(new.dum.vars)){
            #
            #Create name of new dummy variable
            newdummy=paste(list.of.z.factors[i],new.dum.vars[j],sep=".")
            #
            #Assign value of new variable 
            z=cbind(z,assign(newdummy,as.numeric(new.levels.values==new.levels[(j+1)]))) #create new variable based on dummy
            names(z)[dim(z)[2]]=newdummy
          }
          
          #
          #Drop the variable from the dataframe
          z[list.of.z.factors[i]]=NULL
          
          
        }}
      #cat("Transforming all additional covariates variables","\n")
      else{z=scale(z)
      }}}
  z=data.frame(z) 
  names.z=names(z)
  len.z=dim(z)[2]
  row.z=dim(z)[1]
  #cat("The number of observations in covariates are",row.z,"\n")
  if (names(z)[1]=="X1") {
    betanames.z=rep(NA,len.z)
    for (i in 1:len.z){
      betanames.z[i]=paste("z",i,sep="")
      colnames(z)=betanames.z
    }
  } else {
    names.z=names(z)
    betanames.z=names(z)
  }
  
  
  #
  #2.Verify there is no missing data
  if(is.na(sum(y))) stop("Stopping program...missing values in y") 
  if(is.na(sum(sv))) stop("Stopping program...missing values in sv")
  if(is.na(sum(xs))) stop("Stopping program...missing values in xs")
  if(is.na(sum(w))) stop("Stopping program...missing values in w")
  if(is.na(sum(z))) stop("Stopping program...missing values in z")
  
  #
  #Dimensions of data the same
  #Check if other variables have same number of observations as y
  if (!missing(z)){
    if (!all(sapply(list(row.sv, row.xs, row.z, row.w), FUN = identical, row.y)))
    {col.dim=sapply(list(row.sv, row.xs, row.z, row.w), FUN = identical, row.y)
    variable.names=c("Confounding Variables","Instrumental Variables","Covariates","Mismeasured variable")
    for (i in 1:length(variable.names))
    {if (col.dim[i]==FALSE) cat("The",variable.names[i],"has less observations than the response y.")}
    stop("Different dimensions.") 
    }
  } else {
    if (!all(sapply(list(row.sv, row.xs,  row.w), FUN = identical, row.y)))
    {col.dim=sapply(list(row.sv, row.xs,  row.w), FUN = identical, row.y)
    variable.names=c("confounding variables","instrumental variables","mismeasured variable")
    for (i in 1:length(variable.names))
    {if (col.dim[i]==FALSE) cat("The",variable.names[i],"has less observations than the response y.")}
    stop("Different dimensions.") 
    }
    
    
  }
  #
  #Determining number of controls, capm
  case.position=c(which(y==1),length(y))
  if (y[length(y)]==0){
    case.position=c(which(y==1),(length(y)+1))
  }
  length.case.position=length(case.position)
  capm=case.position[2:length.case.position]-case.position[1:(length.case.position-1)]-1
  if (length(unique(capm))!=1)
  {stop("The number of controls is not equal for all stratum.") 
  } else {
    capm=as.numeric(min(capm))
  }

##########FORMAT DATA for analysis
#Need the length of the vectors for the variables to determine
#the number of parameters required
#n=total number of cases
#capm=number of controls for each case

#number of parameters needed for first likelihood
#1 parameter for intercept, parameters for sv and tstar
#1 parameter for misclassification probability if either flag setalpha0.to.0
  #or setalpha1.to.0 is set to true
#otherwise 2 parameters for the misclassification probabilities
  #alpha0 and alpha1
if (setalpha0.to.0==TRUE | setalpha1.to.0==TRUE)
{  
  no.par.one=sum(1,len.sv,len.xs,len.z,1)
} else {
  no.par.one=sum(1,len.sv,len.xs,len.z,2)
}
  cat("The number of gamma and eta parameters is ",no.par.one,". \n",sep="")
#number of parameters needed for the second likelihood
if (!missing(z))
{
  no.par.two=sum(1,len.z)
} else {
  no.par.two=1
}

cat("The number of beta parameters is ",no.par.two,". \n",sep="")
#total number of parameters
n.allpar=no.par.one+no.par.two
cat("The number of total parameters estimating is ",n.allpar,". \n",sep="")


#Used for naive analysis of mismeasured covariate w
#Additionally, naive.estimate.beta is used as a starting point
#for calculating beta
myv=rep(rnorm(n), each=(capm+1))
myst= rep(1:n, each=(capm+1)) 
if(!missing(z)){
  covariates=paste0("z$",names.z,collapse= "+")
  formula = as.formula(paste("Surv(myv, y) ~ ", "w+",covariates,"+strata(myst)"))
} else
{
  formula = as.formula(paste("Surv(myv, y) ~ ", "w+","strata(myst)"))
}

#Used for baseline comparison of results
outw=coxph(formula)
naive.estimate.beta=as.numeric(outw$coef)

#Save results of naive estimates
outw.results=cbind(summary(outw)$coef,summary(outw,conf.int=(1-as.numeric(alpha)))$conf.int[,3:4])

#Used for starting values for gamma and eta in nlglik1.r function
#Additionally, naive.estimate.beta is used as a starting point
#for calculating beta
if(!missing(z)){
  instruments=paste0("xs$",names.xs,collapse= "+")
  confounders=paste0("sv$",names.sv,collapse= "+")
  covariates=paste0("z$",names.z,collapse= "+")
  formula2 = as.formula(paste("w ~ ", instruments,"+",confounders,"+",covariates))
} else
{
  formula2 = as.formula(paste("w ~ ", instruments,"+",confounders))
}
#Used for the starting value in estimating parameters of eta and gamma
startingvalues=glm(formula2,family=binomial)$coef
if (setalpha0.to.0==TRUE | setalpha1.to.0==TRUE)
{  
  startingvalues2=c(startingvalues,0)
} else {
  startingvalues2=c(startingvalues,0,0)
}

#
#
#
##############Begin estimation###############
#Function part1
cat("Conducting estimation of efficient approach...","\n")
if (setalpha0.to.0==TRUE) {
part1.results=.m4noa0(capm=capm, y=y,sv=sv,w=w,xs=xs,z=z,len.xs=len.xs,len.sv=len.sv,len.z=len.z, 
                   no.par.one=no.par.one,no.par.two=no.par.two,startingvalues=startingvalues2,
                   naive.estimate.beta=naive.estimate.beta)  
} else if (setalpha1.to.0==TRUE) {
part1.results=.m4noa1(capm=capm, y=y,sv=sv,w=w,xs=xs,z=z,len.xs=len.xs,len.sv=len.sv,len.z=len.z, 
                   no.par.one=no.par.one,no.par.two=no.par.two,startingvalues=startingvalues2,
                   naive.estimate.beta=naive.estimate.beta)
} else {
part1.results=.m4(capm=capm, y=y,sv=sv,w=w,xs=xs,z=z,len.xs=len.xs,len.sv=len.sv,len.z=len.z, 
                  no.par.one=no.par.one,no.par.two=no.par.two,startingvalues=startingvalues2,
                  naive.estimate.beta=naive.estimate.beta)
}
par.set2=part1.results$par.set.m4
beta.m4=part1.results$beta.m4
var.cov.m4=part1.results$var.cov.m4
par.sd.m4=part1.results$par.set.sd.m4
cat("Done","\n")

#
#
#
##Part 6 - Calculate alpha level, Wald CI, and pvalues for the parameter estimates
alpha=as.numeric(alpha)
cat("Now calculating confidence intervals and p-values with an alpha=",alpha," level of significance.","\n",sep="")
#Create labels for results
param.labels=c("w",betanames.z)

#Calculate results for tables 

cat("Estimation of parameters and standard errors is complete.","\n")
cat("\n")
cat("The results of the analysis of the 1:",capm," matched case control data where the exposure variable is misclassified:","\n",sep="")
cat("\n")

cat("\n")
cat("Analysis using the efficient estimator approach:","\n")
cat("\n")

if (setalpha0.to.0==TRUE) {
final.results2=.finalpartnoa0(par.set=par.set2,beta=beta.m4,no.par.one=no.par.one,
                         par.sd=par.sd.m4,var.cov=var.cov.m4,alpha=alpha,
                         param.labels=param.labels,n.allpar=n.allpar)
} else if (setalpha1.to.0==TRUE) {
final.results2=.finalpartnoa1(par.set=par.set2,beta=beta.m4,no.par.one=no.par.one,
                           par.sd=par.sd.m4,var.cov=var.cov.m4,alpha=alpha,
                           param.labels=param.labels,n.allpar=n.allpar)
  
} else {
final.results2=.finalpart(par.set=par.set2,beta=beta.m4,no.par.one=no.par.one,
                           par.sd=par.sd.m4,var.cov=var.cov.m4,alpha=alpha,
                           param.labels=param.labels,n.allpar=n.allpar)
}

print(final.results2,quote=FALSE,row.names=FALSE)


cat("\n")
cat("Naive analysis using the misclassified variable:","\n")
cat("\n")

print(outw.results,quote=FALSE,row.names=FALSE)


return(list(efficient.results=final.results2,naive.results=outw.results))

#End function
}
