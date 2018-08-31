\name{matcdata}
\alias{matcdata}

\title{
A sample dataset representing matched case control data.  
}

\description{
An example data set with 1 binary response (y), 1 stratification variable (sv), 
1 mismeasured binary exposure (w), the true binary exposure (tr), 1 covariate (z),
and 1 instrumental variable (xs) for use with the meiv function.
}

\usage{
data(matcdata)
}

\seealso{
\code{\link{meiv}}
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
\keyword{ Dataset}
\keyword{ bias}
