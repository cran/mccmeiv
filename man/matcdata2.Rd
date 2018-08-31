\name{matcdata2}
\alias{matcdata2}

\title{
A sample dataset representing matched case control data.  
}

\description{
An example data set with 1 binary response (y), 1 stratification variable (sv), 
1 mismeasured binary exposure (w), the true binary exposure (tr), 
and 1 instrumental variable (xs) with 3 levels (1,2,3) 
for use with the meiv function.

To be used with the 
}

\usage{
data(matcdata2)
}

\seealso{
\code{\link{meiv}}
\code{\link{matcdata}}
\code{\link{matcdatamult}}
}

\examples{
data(matcdata2)
out3=with(matcdata2, meiv(y=y,sv=sv,xs=data.frame(xs=xs),w=w,xs.factor="xs"))
}

\keyword{ Mismeasured}
\keyword{ Matched Cased Control}
\keyword{ Instrumental}
\keyword{ Dataset}
\keyword{ bias}
