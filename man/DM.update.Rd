\name{DM.update}
\alias{DM.update}
\title{FIR.DM updating function}
\usage{
  DM.update(data.train, rule.data.num, miu.rule, func.tsk,
    varinp.mf, step.size = 0.01, def)
}
\arguments{
  \item{data.train}{a matrix (m x n) of normalized data,
  where m is the number of instances and n is the number of
  variables; the last column is the output variable.}

  \item{rule.data.num}{a matrix containing the rulebase.
  Its elements are integers, see \code{\link{rulebase}}.}

  \item{miu.rule}{a matrix with the degrees of rules which
  is a result of the \code{\link{inference}}.}

  \item{func.tsk}{a matrix of parameters of the functions
  on the consequent part of the Takagi Sugeno Kang model.}

  \item{varinp.mf}{a matrix of parameters of the membership
  functions of the input variables.}

  \item{step.size}{the step size of the descent method,
  between 0 and 1.}

  \item{def}{a matrix which is obtained from the
  defuzzification. Please have a look at
  \code{\link{defuzzifier}}.}
}
\description{
  The role of this function is to update the parameters of
  the fuzzy inference rules by descent method (FIR.DM).
  This function is called by the main function of the
  FIR.DM method, \code{\link{FIR.DM}}.
}
\seealso{
  \code{\link{frbs.learn}}, \code{\link{predict}}, and
  \code{\link{FIR.DM}}.
}

