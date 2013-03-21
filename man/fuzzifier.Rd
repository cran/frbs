\name{fuzzifier}
\alias{fuzzifier}
\title{Transform from crisp set into fuzzy terms}
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
  matrix is (5, n) where n is a multiplication the number
  of fuzzy terms/labels and the number of input variables.
  The rows of the matrix represent: The first row is the
  type of membership function, where 1 means
  \code{TRIANGLE}, 2 means \code{TRAPEZOID} in left side, 3
  means \code{TRAPEZOID} in right side, 4 means
  \code{TRAPEZOID} in the middle, 5 means \code{GAUSSIAN},
  6 means \code{SIGMOID}, and 7 means \code{BELL}. And, the
  second up to fifth row indicate the corner points to
  construct the functions. \itemize{ \item \code{TRIANGLE}
  has three parameters (a, b, c), where b is the center
  point of the TRIANGLE, and a and c are the left and right
  points, respectively. \item \code{TRAPEZOID} has four
  parameters (a, b, c, d). \item \code{GAUSSIAN} has two
  parameters (mean and variance). \item \code{SIGMOID} has
  two parameters (gamma and c). \item \code{BELL} has three
  parameters (a, b, c). }

  For example:

  \code{varinp.mf <-
  matrix(c(2,1,3,2,3,0,30,60,0,40,20,50,80,}

  \code{30,80,40,70,100,60,100,0,0,100,0,100), nrow=5,
  byrow=TRUE)}}
}
\value{
  A matrix of the degree of each fuzzy terms based on the
  shape of the membership functions
}
\description{
  Fuzzification refers to the process of transforming a
  crisp set into fuzzy terms.
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

