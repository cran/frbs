\name{fuzzifier}
\alias{fuzzifier}
\title{Transform from crisp set into fuzzy terms}
\usage{
  fuzzifier(data, num.varinput, num.fvalinput, varinp.mf)
}
\arguments{
  \item{data}{a matrix of data containing numerical
  elements.}

  \item{num.varinput}{number of input variables.}

  \item{num.fvalinput}{the number of labels of the input
  variables.}

  \item{varinp.mf}{a matrix containing the parameters to
  form the membership functions. The dimension of the
  matrix is (5, n) where n is the number of fuzzy
  terms/labels and the number of variables. The rows of the
  matrix represent: The first row is the type of membership
  function, where 1 means triangular, 2 means trapezoid 1a
  (left side), 3 means trapezoid 1b (right side), 4 means
  trapezoid 2 (in the middle), 5 means Gaussian, 6 means
  sigmoid, 7 means generalized bell. The second until fifth
  row indicate the critical points to construct the
  functions. \itemize{ \item triangular has three
  parameters (a, b, c), where b is the center point of the
  triangular, and a and c are the left and right points,
  respectively. \item trapezoid has four parameters (a, b,
  c, d). \item Gaussian has two parameters (mean and
  variance). \item sigmoid has two parameters (gamma and
  c). \item generalized bell has three parameters (a, b,
  c). }

  For example:

  \code{varinp.mf <-
  matrix(c(2,1,3,2,3,0,30,60,0,40,20,50,80,}

  \code{30,80,40,70,100,60,100,0,0,100,0,100), nrow=5,
  byrow=TRUE)}}
}
\value{
  A matrix of the degree of each fuzzy term based on the
  shape of the membership functions
}
\description{
  Fuzzification refers to the process of transforming a
  crisp set into fuzzy terms.
}
\details{
  In this function, there are five shapes of membership
  functions implemented, namely triangular, trapezoid,
  Gaussian, sigmoid, and generalized bell.
}
\seealso{
  \code{\link{defuzzifier}}, \code{\link{rulebase}}, and
  \code{\link{inference}}
}

