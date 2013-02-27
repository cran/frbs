\name{HGD.update}
\alias{HGD.update}
\title{FS.HGD updating function}
\usage{
  HGD.update(data.train, miu.rule, func.tsk, varinp.mf,
    step.size = 0.01, def)
}
\arguments{
  \item{data.train}{a matrix(m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables; the last column is the
  output variable.}

  \item{miu.rule}{a matrix with the degrees of rules which
  is the result of the \code{\link{inference}}.}

  \item{func.tsk}{a matrix of parameters of the function on
  the consequent part using the Takagi Sugeno Kang model.
  See \code{\link{rulebase}}.}

  \item{varinp.mf}{a matrix of parameters of membership
  functions of the input variables.}

  \item{step.size}{a real number between 0 and 1
  representing the step size of the gradient descent.}

  \item{def}{a matrix which is obtained by the
  \code{\link{defuzzifier}}.}
}
\description{
  The role of this function is to update parameters within
  the simplified TSK fuzzy rule generation method using
  heuristics and the gradient descent method (FS.HGD). This
  function is called by the main function of the FS.HGD
  method, see \code{\link{FS.HGD}}.
}
\seealso{
  \code{\link{FS.HGD}}
}

