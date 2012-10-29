\name{MSGFS.test}
\alias{MSGFS.test}
\title{MSGFS: The prediction phase}
\usage{
  MSGFS.test(object, newdata)
}
\arguments{
  \item{object}{the \code{\link{frbs-object}}.}

  \item{newdata}{a matrix(m x n) of data for the prediction
  process, where m is the number of instances and n is the
  number of input variables.}
}
\value{
  A matrix of predicted values.
}
\description{
  This function is the internal function of the MSGFS
  method to compute the predicted values.
}

