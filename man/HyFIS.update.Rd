\name{HyFIS.update}
\alias{HyFIS.update}
\title{HyFIS updating function}
\usage{
  HyFIS.update(data.train, def, rule, names.varoutput,
    var.mf, miu.rule, num.labels, MF, step.size = 0.001,
    degree.rule)
}
\arguments{
  \item{data.train}{a matrix(m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables; the last column is the
  output variable.}

  \item{def}{matrix of defuzzification results. See
  \code{\link{defuzzifier}}.}

  \item{rule}{fuzzy IF-THEN rules. See
  \code{\link{rulebase}}.}

  \item{names.varoutput}{a list of names of the output
  variable.}

  \item{var.mf}{a matrix of parameters of the membership
  functions. Please see \code{\link{fuzzifier}}.}

  \item{miu.rule}{a matrix of degree of rules which is a
  result of the \code{\link{inference}}.}

  \item{num.labels}{a matrix(1 x n) whose elements
  represent the number of labels (or fuzzy terms), where n
  is the number of variables.}

  \item{MF}{a matrix of parameters of the membership
  functions which is a result of the
  \code{\link{fuzzifier}}.}

  \item{step.size}{a real number, the step size of the
  gradient descent.}

  \item{degree.rule}{a matrix of degrees of rules. See
  \code{\link{frbs-object}}.}
}
\description{
  This function is called by \code{\link{HyFIS}} to update
  the parameters within the HyFIS method.
}
\seealso{
  \code{\link{HyFIS}}
}

