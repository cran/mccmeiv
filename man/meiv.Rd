\name{meiv}
\alias{meiv}

\title{
Estimating Parameters for a Matched Case Control Design with a Mismeasured Exposure using Instrumental Variables 
}
\description{
Applying the methodology from Manuel, Wang, and Sinha (2018), estimates for the parameters for matched case control data with a mismeasured exposure are calculated
through the use of user supplied instrumental variables.
}
\usage{
meiv(y, sv, xs, w, z = NULL,sv.factor=NULL,xs.factor=NULL,
              z.factor=NULL, alpha=0.05,scale=TRUE)
}

\arguments{
  \item{y}{
A vector of the response variable, representing case (y=1) or control (y=0) for each observation. 
For a 1:M matched case-control dataset with n matched sets, y will be a vector of length N=nx(1+M), 
and it is y=rep( c(1, rep(0, M)), n).

}
  \item{sv}{
A data frame or matrix of instrumental variables used as proxies for the mismeasured variable. 
This data frame should have N number of rows. 
}
  \item{xs}{
A data frame or matrix of instrumental variables used as proxies for the mismeasured variable. 
This data frame should have N number of rows. 
}
  \item{w}{
A vector of the mismeasured binary exposure covariates for each observation. 
The length of this vector must be N.
}
  \item{z}{
A data frame or matrix of additional covariates used to study the association between the response 
and the mismeasured  exposure. If not specified, then only the mismeasured covariate is used. 
This data frame should have N number of rows.
}
  \item{sv.factor}{
Specify whether any stratification variables are factors by using the name of the column(s).
For example if there is a factor variable in the user specified sv labeled "Political Affiliation" and a variable called "SES", then
sv.factor="Political Affiliation" or sv.factor=c("Political Affiliation","SES").
}
  \item{xs.factor}{
Specify whether any instrumental variables are factors by using the name of the column(s).
For example if there is an instrument in the user specified xs labeled "Race of Mother" and a variable called "Race of Father", then
 xs.factor="Race of Mother"  or xs.factor=c("Race of Mother","Race of Father").
}
  \item{z.factor}{
Specify whether any additional variables are factors by using the name of the column(s).
For example, if there is a factor variable in the user specified z labeled Season and a factor variable labeled Ethnicity
z.factor="Season", or z.factor=c("Season","Ethnicity").
}
  \item{alpha}{
Specify the level of significance for calculating the 1-alpha confidence intervals of the estimates. For example, the default alpha=0.05
generates a 95\% confidence interval. 
}
  \item{scale}{
By default, all non factors (those that are not specified as either response, mismeasured covariate, or factors specified in xs.factor/sv.factor/z.factor) are automatically
scaled and centered unless the user specifies FALSE.
}
}
\details{
Estimation of the parameters is done in two steps. In the first step the set of parameters gamma, used to model the probability of true exposure status X -
Pr(X=1|SV,XS,Y=0) and the parameters eta0 and eta1, used to model the misclassification probabilities alpha.0=Pr(W=1|X=0,Y=0)=Pr(W=1|X=0,Y=1)
and alpha.1=Pr(W=0|X=1,Y=0)=Pr(W=0|X=1,Y=1) are estimated. This information is then used in the second step to calculate the Pr(Y=1|SV,XS) or Pr(Y=1|SV,XS,Z) if
the user specified additional covariates Z. In both cases, the solution is found via the optim function, using the "L-BFGS-B" method. For the first step, 
the starting values come from an runif(-0.75,0.75) distribution while the starting values for the second step use the naive beta estimates. In both cases,
a Firth penalty (also known as a Jeffrey's prior) is used to help with convergence. After this the variance covariance matrix for these parameters are estimated. The delta method approach is utilized to estimate the standard errors for the misclassification
probabilities alpha.0 and alpha.1. The time of estimation increases as more covariates (z), instrumental variables (xs) and confounding variables (sv) are added to the model.

In order to specify factors, ensure that the variables containing the factor are correctly labeled - note how the variables for z are specified
in out3, found in the matcdata2 help file.

For an explanation of the data and variables used in the examples, please see the associated data files in the see also section below.

}
\value{

\item{final.results }{Provides estimates for the parameters using the instrumental variable approach and the naive 
analysis. Standard errors, p-values, and a 1-alpha Confidence Interval are also included.
}
\item{par.set }{A list of the estimates from the first analysis. The first value represents the intercept of the logistic model 
Pr(X=1|SV,XS,Y=0), then estimates for each stratification variable specified by sv, estimates for each instrument variable specified 
by XS, followed by the estimates eta.0 and eta.1 which are used to calculate the misclassification probabilities alpha.0 and alpha.1
}
\item{beta }{A list of the estimates from the second analysis. The first value represents the corrected exposure variable estimate, and if 
specified, additional estimates for variables specified in Z. 
}
\item{var.cov }{Variance covariance matrix for the estimates in par.set and beta. 
}
\item{naive.estimates }{Estimates calculated by using the mismeasured variable w, and if specified, additional variables Z.
}

}


\references{
Manuel, CM, Wang, S, and Sinha, S (2018). Matched case-control data with a misclassified exposure: 
What can be done with instrumental variables? 
}
\author{
Chris M. Manuel and Samiran Sinha 
}

\seealso{
\code{\link{matcdata}}
\code{\link{matcdata2}}
\code{\link{matcdatamult}}
}

\examples{
data(matcdata)
out1=with(matcdata,meiv(y=y,sv=sv,xs=xs,w=w,z=z,alpha=0.1))
}

\keyword{ Mismeasured}
\keyword{ Matched Cased Control}
\keyword{ Instrumental}
\keyword{ bias}
