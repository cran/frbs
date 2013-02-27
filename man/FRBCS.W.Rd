\name{FRBCS.W}
\alias{FRBCS.W}
\title{FRBCS.W model building}
\usage{
  FRBCS.W(range.data, data.train, label.inp, num.class,
    type.mf)
}
\arguments{
  \item{range.data}{a matrix(2 x n) containing the range of
  the normalized data, where n is the number of variables,
  and first and second rows are the minimum and maximum
  values, respectively.}

  \item{data.train}{a matrix(m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables; the last column is the
  output variable.}

  \item{label.inp}{a matrix(1 x n) whose elements represent
  the number of labels (fuzzy terms), where n is the number
  of variables.}

  \item{num.class}{an integer number representing the
  number of labels (fuzzy terms).}

  \item{type.mf}{the type of the shape of the membership
  functions.}
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

