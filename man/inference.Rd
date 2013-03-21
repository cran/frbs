\name{inference}
\alias{inference}
\title{The process of fuzzy reasoning}
\usage{
  inference(MF, rule, names.varinput, type.tnorm,
    type.snorm)
}
\arguments{
  \item{MF}{a matrix of the degrees of membership functions
  which is a result of the \code{\link{fuzzifier}}.}

  \item{rule}{a matrix or list of fuzzy IF-THEN rules. See
  \code{\link{rulebase}}.}

  \item{names.varinput}{a list of names of the input
  variables.}

  \item{type.tnorm}{a value which represents the type of
  t-norm to be used: \itemize{ \item \code{1} or \code{MIN}
  means standard t-norm: min(x1, x2). \item \code{2} or
  \code{HAMACHER} means Hamacher product: (x1 * x2)/(x1 +
  x2 - x1 * x2). \item \code{3} or \code{YAGER} means Yager
  class (with tao = 1): 1- min(1, ((1 - x1) + (1 - x2))).
  \item \code{4} or \code{PRODUCT} means product: (x1 *
  x2). \item \code{5} or \code{BOUNDED} means bounded
  product: max(0, x1 + x2 - 1). }}

  \item{type.snorm}{a value which represents the type of
  s-norm to be used: \itemize{ \item \code{1} or \code{MAX}
  means standard s-norm: max(x1, x2). \item \code{2} or
  \code{HAMACHER} means Hamacher sum: (x1 + x2 - 2x1 * x2)
  / 1 - x1 * x2. \item \code{3} or \code{YAGER} means Yager
  class (with tao = 1): min(1, (x1 + x2)). \item \code{4}
  or \code{SUM} means sum: (x1 + x2 - x1 * x2). \item
  \code{5} or \code{BOUNDED} means bounded sum: min(1, x1 +
  x2). }}
}
\value{
  a matrix of the degrees of the rules.
}
\description{
  Inference refers to the process of fuzzy reasoning.
}
\details{
  In this function, fuzzy reasoning is conducted based on
  Mamdani and Takagi Sugeno Kang model. Furthermore, there
  are some formula for conjunction and disjunction
  operators.

  \bold{The Mamdani model:} A fuzzy system with, e.g., two
  inputs x1 and x2 (antecedents) and a single output y
  (consequent) is described by the following fuzzy IF-THEN
  rule:

  \code{IF x1 is A1 and x2 is A2 THEN y is B}

  where A1 and A2 are the fuzzy sets representing the
  antecent pairs and B is the fuzzy set representing the
  consequent.

  \bold{The Takagi Sugeno Kang model:} Suppose we have two
  inputs x1 and x2 and output y, then the fuzzy IF-THEN
  rule is as follows:

  \code{IF x1 is A1 and x2 is A2 THEN y is y = f(x1, x2)}

  where y = f(x1, x2) is a crisp function in the consequent
  part which is usually a polynomial function, and A1 and
  A2 are the fuzzy sets representing the antecent pairs.

  Futhermore, this function has the following capabilities:
  \itemize{ \item It supports unary operators (not) and
  binary operators (AND and OR). \item It provides
  linguistic hedge (extremely, very, somewhat, slightly).
  \item there are several methods for the t-norm and
  s-norm. }
}
\seealso{
  \code{\link{defuzzifier}}, \code{\link{rulebase}}, and
  \code{\link{fuzzifier}}.
}

