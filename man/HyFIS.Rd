\name{HyFIS}
\alias{HyFIS}
\title{HyFIS model building}
\usage{
  HyFIS(range.data, data.train, num.labels, max.iter = 100,
    range.data.ori, step.size = 0.01)
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

  \item{num.labels}{a matrix(1 x n), whose elements
  represent the number of labels (fuzzy terms); n is the
  number of variables.}

  \item{max.iter}{the maximal number of iterations.}

  \item{range.data.ori}{a matrix containing the ranges of
  the original data.}

  \item{step.size}{step size of the gradient descent
  method.}
}
\description{
  This is the internal function that implements the hybrid
  neural fuzzy inference system (HyFIS). Users do not need
  to call it directly, but just use
  \code{\link{frbs.learn}} and \code{\link{predict}}
}
\details{
  This method was proposed by J. Kim and N. Kasabov. There
  are two phases in this method for learning, namely the
  knowledge acquisition module and the structure and
  parameter learning. The knowledge acquition module uses
  the techniques of Wang and Mendel. The learning of
  structure and parameters is a supervised learning method
  using gradient descent-based learning algorithms. This
  function generates a model which consists of a rule
  database and parameters of the membership functions. The
  rules of HyFIS use the Mamdani model on the antecedent
  and consequent parts. Futhermore, HyFIS uses a Gaussian
  membership function. So, there are two kinds of
  parameters that are optimized, mean and variance of the
  Gaussian function.
}
\references{
  J. Kim and N. Kasabov, "HyFIS: adaptive neuro-fuzzy
  inference systems and their application to nonlinear
  dynamical systems", Neural Networks 12, 1301 - 1319
  (1999).
}
\seealso{
  \code{\link{HyFIS.update}}, \code{\link{frbs.learn}}, and
  \code{\link{predict}}.
}

