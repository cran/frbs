\name{GFS.Thrift.test}
\alias{GFS.Thrift.test}
\title{GFS.Thrift: The prediction phase}
\usage{
  GFS.Thrift.test(object, newdata)
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
  This function is the internal function of the GFS.Thrift
  method to compute the predicted values.
}

