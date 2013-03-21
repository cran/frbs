\name{GFS.LT.RS.test}
\alias{GFS.LT.RS.test}
\title{GFS.LT.RS: The prediction phase}
\usage{
  GFS.LT.RS.test(object, newdata)
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
  This function is the internal function of the GFS.LT.RS
  method to compute the predicted values.
}

