\name{HGD}
\alias{HGD}
\title{HGD model building}
\usage{
  HGD(range.data, data.train, num.labels, max.iter,
    step.size, alpha.heuristic = 1)
}
\arguments{
  \item{data.train}{a matrix(m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables; the last column is the
  output variable.}

  \item{range.data}{a matrix(2 x n) containing the range of
  the normalized data, where n is the number of variables,
  and first and second rows are the minimum and maximum
  value, respectively.}

  \item{num.labels}{a matrix(1 x n), whose elements
  represent the number of labels (fuzzy terms); n is the
  number of variables.}

  \item{max.iter}{maximal number of iterations.}

  \item{step.size}{step size of the descent method.}

  \item{alpha.heuristic}{a positive real number which is
  the heuristic parameter.}
}
\description{
  This is the internal function that implements the fuzzy
  system using heuristics and gradient descent method
  (HGD). Users do not need to call it directly, but just
  use \code{\link{frbs.learn}} and \code{\link{predict}}.
}
\details{
  This method was proposed by Ken Nozaki, H. Ishibuchi, and
  Hideo Tanaka. It uses fuzzy IF-THEN rules with nonfuzzy
  singletons (i.e. real numbers) in the consequent parts.
  The techniques of Wang and Mendel are implemented to
  generate the antecedent part, while the initial
  consequent part of each rule is determined by the
  weighted mean value of the given training data. Then, the
  gradient descent method updates the value of the
  consequent part. Futhermore, the heuristic value given by
  the user affects the value of weight of each data.
}
\references{
  H. Ichihashi and T. Watanabe, "Learning control system by
  a simplified fuzzy reasoning model," Proc. IPMU'90 417 -
  419 (1990).

  H. Ishibuchi, K. Nozaki, H. Tanaka, Y. Hosaka and M.
  Matsuda, "Empirical study on learning in fuzzy systems by
  rice taste analysis," Fuzzy Set and Systems 64, 129 144
  (1994).
}
\seealso{
  \code{\link{frbs.learn}}, \code{\link{predict}}, and
  \code{\link{HGD.update}}
}

