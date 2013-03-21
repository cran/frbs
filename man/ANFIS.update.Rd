\name{ANFIS.update}
\alias{ANFIS.update}
\title{ANFIS updating function}
\usage{
  ANFIS.update(data.train, range.input, range.output,
    rule.data.num, miu.rule, func.tsk, varinp.mf,
    step.size = 0.01)
}
\arguments{
  \item{data.train}{a matrix(m x n) of normalized data for
  the training process, where m is the number of instances
  and n is the number of variables; the last column is the
  output variable.}

  \item{range.input}{the range of the input variables, as a
  matrix(2 x n).}

  \item{range.output}{the range of the output variable, as
  a matrix(2 x n).}

  \item{rule.data.num}{a matrix containing the rule base in
  integer form.}

  \item{miu.rule}{a matrix with the degrees of rules. See
  \code{\link{inference}}.}

  \item{func.tsk}{a matrix of parameters of the function on
  the consequent part using the Takagi Sugeno Kang model.}

  \item{varinp.mf}{a matrix of parameters of membership
  functions of the input variables.}

  \item{step.size}{a real number between 0 and 1
  representing the step size of the gradient descent.}
}
\description{
  The role of this function is to update parameters in the
  ANFIS method. This function is called by the main
  function of the ANFIS method, \code{\link{ANFIS}}.
}
\seealso{
  \code{\link{ANFIS}}, \code{\link{rulebase}},
  \code{\link{inference}}, \code{\link{fuzzifier}} and
  \code{\link{defuzzifier}}
}

