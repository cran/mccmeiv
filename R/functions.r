#Begin file

#Negative row sum
#This function sums the negative for each multiplicative term of parameter and column
#Like rowSum but does not sum all columns at once
#Mainly for when there are more than one instrument or stratification variable
.nrSum=function(par,columns){
  obsum=0
  ncol=dim(columns)[2]
  npar=length(par)
  if(ncol!=npar){
    print("number of instruments not equal to number of parameters")
  }
  for (i in 1:npar ) {obsum=obsum-par[i]*columns[,i]}
  obsum=round(obsum,digits=7) 
  return(obsum)
} 

#Positive row sum
#This function sums the negative for each multiplicative term of parameter and column
#Like rowSum but does not sum all columns at once
#Mainly for when there are more than one instrument or stratification variable
.prSum=function(par,columns){
  obsum=0
  ncol=dim(columns)[2]
  npar=length(par)
  if(ncol!=npar){
    print("number of instruments not equal to number of parameters")
  }
  for (i in 1:npar ) {obsum=obsum+par[i]*columns[,i]}
  obsum=round(obsum,digits=7) 
  return(obsum)
} 
