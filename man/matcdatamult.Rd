\name{matcdatamult}
\alias{matcdatamult}

\title{
A sample dataset representing matched case control data. Similar to matcdata except with multiple 
stratification/instrumental/miscellaneous variables.
}

\description{
An example data set with 1 binary response (y), 2 stratification variable (sv1 and sv2), 
1 mismeasured binary exposure (w),  2 prognostic factors (z1 and z2),
and 2 instruments (xs1 and xs2) for use with the meiv function.
}

\usage{
data(matcdatamult)
}

\seealso{
\code{\link{meivm3}}
\code{\link{meivm4}}
\code{\link{matcdata}}
}

\examples{
data(matcdatamult)
out=with(matcdatamult,meivm3(y=y, sv=sv1,xs=xs1, w=w,z=cbind(z1,z2)))
}

\keyword{ Mismeasured}
\keyword{ Matched Cased Control}
\keyword{ Instrumental}
\keyword{ Dataset}
\keyword{ bias}
