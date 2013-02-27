\name{defuzzifier}
\alias{defuzzifier}
\title{Defuzzifier to transform from fuzzy terms to crisp values}
\usage{
  defuzzifier(data, rule, range.output,
    names.varoutput = NULL, varout.mf = NULL, miu.rule,
    type.defuz = 1, type.model = 1, func.tsk = NULL)
}
\arguments{
  \item{data}{a matrix(m x n) of data, where m is the
  number of instances and n is the number of variables.}

  \item{rule}{a list or matrix of fuzzy IF-THEN rules, as
  discussed in \code{\link{rulebase}}.}

  \item{range.output}{a matrix(2 x n) containing the range
  of the output data.}

  \item{names.varoutput}{a list for giving names to the
  fuzzy terms. See \code{\link{rulebase}}.}

  \item{varout.mf}{a matrix constructing the membership
  function of the output variable.  See
  \code{\link{fuzzifier}}.}

  \item{miu.rule}{the results of the inference module. See
  \code{\link{inference}}.}

  \item{type.defuz}{the type of defuzzification to be used,
  where 1 means weighted average method, and 2, 3, 4 and 5
  mean first, last, mean maxima and modified COG,
  respectively.}

  \item{type.model}{the type of the model that will be used
  in the simulation.  Here, 1 or 2 means we use Mamdani or
  Takagi Sugeno Kang model, respectively.}

  \item{func.tsk}{a matrix used to build the linear
  equation for the consequent part if we are using Takagi
  Sugeno Kang. See also \code{\link{rulebase}}.}
}
\value{
  A matrix of crisp values
}
\description{
  Defuzzification is a transformation that extracts the
  crisp values from the fuzzy terms.
}
\details{
  In this function, there exist two kinds of models which
  are based on the Mamdani model and the Takagi Sugeno Kang
  model. For the Mamdani model there are five methods for
  defuzzifying a fuzzy term A of a universe of discourse Z.
  They are as follows: \enumerate{ \item weighted average
  method \item first of maxima \item last of maxima \item
  mean of maxima \item modified COG }
}
\seealso{
  \code{\link{fuzzifier}}, \code{\link{rulebase}}, and
  \code{\link{inference}}
}

