\name{defuzzifier}
\alias{defuzzifier}
\title{Defuzzifier to transform from linguistic terms to crisp values}
\usage{
  defuzzifier(data, rule = NULL, range.output = NULL,
    names.varoutput = NULL, varout.mf = NULL, miu.rule,
    type.defuz = NULL, type.model = "TSK", func.tsk = NULL)
}
\arguments{
  \item{data}{a matrix (\eqn{m \times n}) of data, where
  \eqn{m} is the number of instances and \eqn{n} is the
  number of variables.}

  \item{rule}{a list or matrix of fuzzy IF-THEN rules, as
  discussed in \code{\link{rulebase}}.}

  \item{range.output}{a matrix (\eqn{2 \times n})
  containing the range of the output data.}

  \item{names.varoutput}{a list for giving names to the
  linguistic terms. See \code{\link{rulebase}}.}

  \item{varout.mf}{a matrix constructing the membership
  function of the output variable.  See
  \code{\link{fuzzifier}}.}

  \item{miu.rule}{the results of the inference module. See
  \code{\link{inference}}.}

  \item{type.defuz}{the type of defuzzification to be used
  as follows.  \itemize{ \item \code{1} or \code{WAM} means
  weighted average method, \item \code{2} or
  \code{FIRST.MAX} means first maxima, \item \code{3} or
  \code{LAST.MAX} means last maxima, \item \code{4} or
  \code{MEAN.MAX} means mean maxima, \item \code{5} or
  \code{COG} means modified center of gravity (COG).  }}

  \item{type.model}{the type of the model that will be used
  in the simulation.  Here, \code{1} or \code{MAMDANI} and
  \code{2} or \code{TSK} means we use Mamdani or Takagi
  Sugeno Kang model, respectively.}

  \item{func.tsk}{a matrix used to build the linear
  equation for the consequent part if we are using Takagi
  Sugeno Kang. See also \code{\link{rulebase}}.}
}
\value{
  A matrix of crisp values
}
\description{
  Defuzzification is a transformation that extracts the
  crisp values from the linguistic terms.
}
\details{
  In this function, there exist two kinds of models which
  are based on the Mamdani and Takagi Sugeno Kang model.
  For the Mamdani model there are five methods for
  defuzzifying a linguistic term \eqn{A} of a universe of
  discourse \eqn{Z}. They are as follows: \enumerate{ \item
  weighted average method (\code{WAM}). \item first of
  maxima (\code{FIRST.MAX}). \item last of maxima
  (\code{LAST.MAX}) \item mean of maxima (\code{MEAN.MAX}).
  \item modified center of gravity (\code{COG}). }
}
\seealso{
  \code{\link{fuzzifier}}, \code{\link{rulebase}}, and
  \code{\link{inference}}
}

