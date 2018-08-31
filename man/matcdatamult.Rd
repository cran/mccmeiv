\name{matcdatamult}
\alias{matcdatamult}

\title{
A sample dataset representing matched case control data. Similar to matcdata except with multiple 
stratification/instrumental/miscellaneous variables.
}

\description{
An example data set with 1 binary response (y), 2 stratification variable (sv1 and sv2), 
1 mismeasured binary exposure (w), the true binary exposure (tr), 2 covariates (z1 and z2),
and 2 instruments (xs1 and xs2) for use with the meiv function.
}

\usage{
data(matcdatamult)
}

\seealso{
\code{\link{meiv}}
\code{\link{matcdata}}
\code{\link{matcdata2}}
}

\examples{
data(matcdatamult)
out2=with(matcdatamult,meiv(y=y, sv=cbind(sv1,sv2),xs=xs1, w=w))
}

\keyword{ Mismeasured}
\keyword{ Matched Cased Control}
\keyword{ Instrumental}
\keyword{ Dataset}
\keyword{ bias}
