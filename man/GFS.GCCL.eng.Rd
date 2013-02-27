\name{GFS.GCCL.eng}
\alias{GFS.GCCL.eng}
\title{GFS.GCCL.test: The prediction phase}
\usage{
  GFS.GCCL.eng(object, newdata)
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
  This function is the internal function of the GFS.GCCL
  and FH.GBML method to compute the predicted values.
}

