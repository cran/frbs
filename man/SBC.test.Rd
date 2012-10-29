\name{SBC.test}
\alias{SBC.test}
\title{SBC prediction phase}
\usage{
  SBC.test(object, newdata)
}
\arguments{
  \item{object}{the \code{\link{frbs-object}}.}

  \item{newdata}{a matrix(m x n) of data for the prediction
  process, where m is the number of instances and n is the
  number of input variables.}
}
\description{
  This function is the internal function of the SBC method
  to compute the predicted values.
}
\seealso{
  \code{\link{SBC}}
}

