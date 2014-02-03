\name{FRBCS.W}
\alias{FRBCS.W}
\title{FRBCS.W model building}
\usage{
  FRBCS.W(range.data, data.train, num.labels, num.class,
    type.mf, type.tnorm = "MIN", type.snorm = "MAX",
    type.implication.func = "ZADEH")
}
\arguments{
  \item{range.data}{a matrix (\eqn{2 \times n}) containing
  the range of the normalized data, where \eqn{n} is the
  number of variables, and first and second rows are the
  minimum and maximum values, respectively.}

  \item{data.train}{a matrix (\eqn{m \times n}) of
  normalized data for the training process, where \eqn{m}
  is the number of instances and \eqn{n} is the number of
  variables; the last column is the output variable. Note
  the data must be normalized between 0 and 1.}

  \item{num.labels}{a matrix (\eqn{1 \times n}), whose
  elements represent the number of labels (linguistic
  terms); \eqn{n} is the number of variables.}

  \item{num.class}{an integer number representing the
  number of labels (linguistic terms).}

  \item{type.mf}{the type of the shape of the membership
  functions.}

  \item{type.tnorm}{the type of t-norm. See
  \code{\link{inference}}.}

  \item{type.snorm}{the type of s-norm. See
  \code{\link{inference}}.}

  \item{type.implication.func}{the type of implication
  function. See \code{\link{WM}}.}
}
\description{
  This is the internal function that implements the fuzzy
  rule-based classification system with weight factor
  (FRBCS.W). It is used to solve classification tasks.
  Users do not need to call it directly, but just use
  \code{\link{frbs.learn}} and \code{\link{predict}}. This
  method is suitable only for classification problems.
}
\details{
  This method is adopted from Hisao Ishibuchi and Tomoharu
  Nakashima's paper. Each fuzzy IF-THEN rule consists of
  antecedent linguistic values and a single consequent
  class with certainty grades (weights). The antecedent
  part is determined by a grid-type fuzzy partition from
  the training data. The consequent class is defined as the
  dominant class in the fuzzy subspace corresponding to the
  antecedent part of each fuzzy IF-THEN rule and the
  certainty grade is calculated from the ratio among the
  consequent class. A class of the new instance is
  determined by the consequent class of the rule with the
  maximal product of the compatibility grade and the
  certainty grade.
}
\references{
  H. Ishibuchi and T. Nakashima, "Effect of rule weights in
  fuzzy rule-based classification systems", IEEE
  Transactions on Fuzzy Systems, vol. 1, pp. 59 - 64
  (2001).
}
\seealso{
  \code{\link{FRBCS.eng}}, \code{\link{frbs.learn}}, and
  \code{\link{predict}}
}

