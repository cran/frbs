\name{frbs.eng}
\alias{frbs.eng}
\title{The prediction phase}
\usage{
  frbs.eng(object, newdata)
}
\arguments{
  \item{object}{the \code{\link{frbs-object}}.}

  \item{newdata}{a matrix (\eqn{m \times n}) of data for
  the prediction process, where \eqn{m} is the number of
  instances and \eqn{n} is the number of input variables.}
}
\value{
  A list with the following items: \item{\code{rule}}{the
  fuzzy IF-THEN rules} \item{\code{varinp.mf}}{a matrix to
  generate the shapes of the membership functions for the
  input variables} \item{\code{MF}}{a matrix of the degrees
  of the membership functions} \item{\code{miu.rule}}{a
  matrix of the degrees of the rules}
  \item{\code{func.tsk}}{a matrix of the Takagi Sugeno Kang
  model for the consequent part of the fuzzy IF-THEN rules}
  \item{\code{predicted.val}}{a matrix of the predicted
  values}
}
\description{
  This function is one of the main internal functions of
  the package. It determines the values within the
  prediction phase.
}
\details{
  This function involves four different processing steps on
  fuzzy rule-based systems. Firstly, the rulebase (see
  \code{\link{rulebase}}) validates the consistency of the
  fuzzy IF-THEN rules form. Then, the fuzzification (see
  \code{\link{fuzzifier}}) transforms crisp values into
  linguistic terms. Next, the inference calculates the
  degree of rule strengths using the t-norm and the s-norm.
  Finally, the defuzzification process calculates the
  results of the model using the Mamdani or the Takagi
  Sugeno Kang model.
}
\seealso{
  \code{\link{fuzzifier}}, \code{\link{rulebase}},
  \code{\link{inference}} and \code{\link{defuzzifier}}.
}

