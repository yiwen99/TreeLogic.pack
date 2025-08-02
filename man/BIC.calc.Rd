\name{BIC.calc}
\alias{BIC.calc}
\title{BIC.calc}
\usage{
BIC.calc(y, pred)
}
\arguments{
  \item{y}{a response vector}
  \item{pred}{the design matrix}
}
\description{
This function is a helper function that calculates BIC. It is used in the Bottom_up_selection and Top_down selection function.
}
\examples{
y <- rnorm(n = 10, mean = 0, sd = 1)
pred <- cbind(rnorm(n = 10, mean = -1, sd = 2), rnorm(n = 10, mean = 1, sd = 1))
BIC.calc(y, pred)
}
\author{
Yiwen Chen
}
