\name{DENFIS.eng}
\alias{DENFIS.eng}
\title{DENFIS prediction function}
\usage{
  DENFIS.eng(object, newdata)
}
\arguments{
  \item{object}{the frbs model. See
  \code{\link{frbs-object}}.}

  \item{newdata}{a matrix(m x n) of data for the prediction
  process, where m is the number of instances and n is the
  number of input variables.}
}
\value{
  a matrix of predicted values
}
\description{
  This function is an internal function for the prediction
  phase using the DENFIS method. The user should use this
  function not directly, but with calling
  \code{\link{predict}}.
}
\seealso{
  \code{\link{DENFIS}}
}

