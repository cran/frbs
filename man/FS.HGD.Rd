\name{FS.HGD}
\alias{FS.HGD}
\title{FS.HGD model building}
\usage{
  FS.HGD(data.train, num.labels, max.iter = 100,
    step.size = 0.01, alpha.heuristic = 1,
    type.tnorm = "MIN", type.snorm = "MAX",
    type.implication.func = "ZADEH")
}
\arguments{
  \item{data.train}{a matrix (m x n) of normalized data for
  the training process, where m is the number of instances
  and n is the number of variables; the last column is the
  output variable.  Note the data must be normalized
  between 0 and 1.}

  \item{num.labels}{a matrix (1 x n), whose elements
  represent the number of labels (fuzzy terms); n is the
  number of variables.}

  \item{max.iter}{maximal number of iterations.}

  \item{step.size}{step size of the descent method.}

  \item{alpha.heuristic}{a positive real number which is
  the heuristic parameter.}

  \item{type.tnorm}{the type of t-norm. For more detail,
  please have a look at \code{\link{inference}}.}

  \item{type.snorm}{the type of s-norm. For more detail,
  please have a look at \code{\link{inference}}.}

  \item{type.implication.func}{a value representing type of
  implication function. For more detail, please have a look
  at \code{\link{WM}}.}
}
\description{
  This is the internal function that implements the
  simplified TSK fuzzy rule generation method using
  heuristics and gradient descent method (FS.HGD). It is
  used to solve regression tasks. Users do not need to call
  it directly, but just use \code{\link{frbs.learn}} and
  \code{\link{predict}}.
}
\details{
  This method was proposed by Ken Nozaki, H. Ishibuchi, and
  Hideo Tanaka. It uses fuzzy IF-THEN rules with nonfuzzy
  singletons (i.e. real numbers) in the consequent parts.
  The techniques of space partition are implemented to
  generate the antecedent part, while the initial
  consequent part of each rule is determined by the
  weighted mean value of the given training data. Then, the
  gradient descent method updates the value of the
  consequent part. Futhermore, the heuristic value given by
  the user affects the value of weight of each data.
}
\references{
  H. Ishibuchi, K. Nozaki, H. Tanaka, Y. Hosaka and M.
  Matsuda, "Empirical study on learning in fuzzy systems by
  rice taste analysis", Fuzzy Set and Systems, vol. 64, no.
  2, pp. 129 - 144 (1994).
}
\seealso{
  \code{\link{frbs.learn}}, \code{\link{predict}}, and
  \code{\link{HGD.update}}
}

