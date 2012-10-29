\name{DM}
\alias{DM}
\title{DM model building}
\usage{
  DM(range.data, data.train, num.labels, max.iter,
    step.size)
}
\arguments{
  \item{range.data}{a matrix(2 x n) containing the range of
  the normalized data, where n is the number of variables,
  and first and second rows are the minimum and maximum
  value, respectively.}

  \item{data.train}{a matrix(m x n) of data for training,
  where m is the number of instances and n is the number of
  variables. The last column is the output variable.}

  \item{num.labels}{a matrix(1 x n) whose elements
  represent the number of labels (fuzzy terms), where n is
  the number of variables.}

  \item{max.iter}{the maximal number of iterations.}

  \item{step.size}{the step size of the descent method,
  between 0 and 1.}
}
\description{
  This is the internal function that implements the fuzzy
  inference rules by descent method (DM). Users do not need
  to call it directly, but just use
  \code{\link{frbs.learn}} and \code{\link{predict}}.
}
\details{
  This method was proposed by Hiroyoshi Nomura, Isao
  Hayashi, and Noboru Wakami. DM uses simplified fuzzy
  reasoning where the consequent part is a real number (a
  particular case within the Takagi Sugeno Kang model),
  while the membership function on the antecedent part is
  expressed by an isosceles triangle. So, in the learning
  phase, DM updates three parameters which are center and
  width of the triangular and a real number on the
  consequent part using a descent method.
}
\references{
  H. Nomura, I. Hayashi and N. Wakami, "A learning method
  of fuzzy inference rules by descent method," Proc.
  FUZZ-IEEE'92 pp. 203 - 210 (1992).
}
\seealso{
  \code{\link{DM.update}}
}

