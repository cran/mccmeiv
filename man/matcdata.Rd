\name{matcdata}
\alias{matcdata}

\title{
A sample dataset representing matched case control data.  
}

\description{
An example data set with 1 binary response (y), 1 stratification variable (sv), 
1 mismeasured binary exposure (w), 1 prognostic factor (z),
and 1 instrumental variable (xs) for use with the meivm3 or meivm4 functions.
}

\usage{
data(matcdata)
}

\seealso{
\code{\link{meivm3}}
\code{\link{meivm4}}
\code{\link{matcdatamult}}
}

\examples{
data(matcdata)
out1=with(matcdata,meivm4(y=y,sv=sv,xs=xs,w=w,z=z,alpha=0.1))
}

\keyword{ Mismeasured}
\keyword{ Matched Cased Control}
\keyword{ Instrumental}
\keyword{ Dataset}
\keyword{ bias}
