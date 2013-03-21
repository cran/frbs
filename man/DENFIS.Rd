\name{DENFIS}
\alias{DENFIS}
\title{DENFIS model building}
\usage{
  DENFIS(data.train, range.data.ori, Dthr = 0.1,
    max.iter = 100, step.size = 0.01, d = 2)
}
\arguments{
  \item{data.train}{a matrix (m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables (input and output
  variables).}

  \item{range.data.ori}{a matrix (2 x n) containing the
  range of the data, where n is the number of variables,
  and first and second rows are the minimum and maximum
  values, respectively.}

  \item{Dthr}{the threshold value for the evolving
  clustering method (ECM), between 0 and 1.}

  \item{max.iter}{the maximal number of iterations.}

  \item{step.size}{the step size of the least squares
  method, between 0 and 1.}

  \item{d}{a parameter for the width of the triangular
  membership function.}
}
\description{
  This is the internal function that implements the dynamic
  evolving neural-fuzzy inference system (DENFIS). It is
  used to handle regression tasks. Users do not need to
  call it directly, but just use \code{\link{frbs.learn}}
  and \code{\link{predict}}.
}
\details{
  This method was proposed by Nikola K. Kasabov and Qun
  Song. There are several steps in this method that are to
  determine the cluster centers using the evolving
  clustering method (ECM), to partition the input space and
  to find optimal parameters on the consequent part (Takagi
  Sugeno Kang model) for the IF-THEN rule using a least
  squares estimator.

  ECM is a distance-based clustering method which is
  determined by a threshold value, \code{Dthr}. This
  parameter influences how many clusters are created. In
  the beginning of the clustering process, the first
  instance from the training data is chosen to be a cluster
  center, and the determining radius is set to zero.
  Afterwards, using the next instance, cluster centers and
  radius are changed based on certain mechanisms of ECM
  (please see \code{\link{ECM}}). All of the cluster
  centers are then obtained after evaluating all the
  training data. The next step is to update the parameters
  on the consequent part with the assumption that the
  antecedent part which we got from ECM is fixed. Actually,
  ECM can perform well as an online clustering method, but
  in this package it is used in an offline mode.
}
\references{
  N. K. Kasabov and Q. Song, "DENFIS: dynamic evolving
  neural-fuzzy inference system and its Application for
  time-series prediction", IEEE Transactions on Fuzzy
  Systems, vol. 10, no. 2, pp. 144 - 154 (2002).
}
\seealso{
  \code{\link{DENFIS.eng}}, \code{\link{frbs.learn}}, and
  \code{\link{predict}}
}

