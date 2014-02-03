\name{fuzzifier}
\alias{fuzzifier}
\title{Transform from crisp set into linguistic terms}
\usage{
  fuzzifier(data, num.varinput, num.labels.input,
    varinp.mf)
}
\arguments{
  \item{data}{a matrix of data containing numerical
  elements.}

  \item{num.varinput}{number of input variables.}

  \item{num.labels.input}{the number of labels of the input
  variables.}

  \item{varinp.mf}{a matrix containing the parameters to
  form the membership functions. The dimension of the
  matrix is (\eqn{5, n}) where \eqn{n} is a multiplication
  the number of linguistic terms/labels and the number of
  input variables. The rows of the matrix represent: The
  first row is the type of membership function, where 1
  means \code{TRIANGLE}, 2 means \code{TRAPEZOID} in left
  side, 3 means \code{TRAPEZOID} in right side, 4 means
  \code{TRAPEZOID} in the middle, 5 means \code{GAUSSIAN},
  6 means \code{SIGMOID}, and 7 means \code{BELL}. And, the
  second up to fifth row indicate the corner points to
  construct the functions. \itemize{ \item \code{TRIANGLE}
  has three parameters (\eqn{a, b, c}), where \eqn{b} is
  the center point of the \code{TRIANGLE}, and \eqn{a} and
  \eqn{c} are the left and right points, respectively.
  \item \code{TRAPEZOID} has four parameters (\eqn{a, b, c,
  d}). \item \code{GAUSSIAN} has two parameters (\eqn{mean}
  and \eqn{variance}). \item \code{SIGMOID} has two
  parameters (\eqn{\gamma} and \eqn{c}). \item \code{BELL}
  has three parameters (\eqn{a, b, c}). }

  For example:

  \code{varinp.mf <-
  matrix(c(2,1,3,2,3,0,30,60,0,40,20,50,80,}

  \code{30,80,40,70,100,60,100,0,0,100,0,100), nrow=5,
  byrow=TRUE)}}
}
\value{
  A matrix of the degree of each linguistic terms based on
  the shape of the membership functions
}
\description{
  Fuzzification refers to the process of transforming a
  crisp set into linguistic terms.
}
\details{
  In this function, there are five shapes of membership
  functions implemented, namely \code{TRIANGLE},
  \code{TRAPEZOID}, \code{GAUSSIAN}, \code{SIGMOID}, and
  \code{BELL}.
}
\seealso{
  \code{\link{defuzzifier}}, \code{\link{rulebase}}, and
  \code{\link{inference}}
}

