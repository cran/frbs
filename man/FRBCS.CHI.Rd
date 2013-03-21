\name{FRBCS.CHI}
\alias{FRBCS.CHI}
\title{FRBCS.CHI model building}
\usage{
  FRBCS.CHI(range.data, data.train, num.labels, num.class,
    type.mf = "TRIANGLE", type.tnorm = "MIN",
    type.snorm = "MAX", type.implication.func = "ZADEH")
}
\arguments{
  \item{range.data}{a matrix (2 x n) containing the range
  of the normalized data, where n is the number of
  variables, and first and second rows are the minimum and
  maximum values, respectively.}

  \item{data.train}{a matrix (m x n) of normalized data for
  the training process, where m is the number of instances
  and n is the number of variables; the last column is the
  output variable. Note the data must be normalized between
  0 and 1.}

  \item{num.labels}{a matrix (1 x n), whose elements
  represent the number of labels (fuzzy terms); n is the
  number of variables.}

  \item{num.class}{an integer number representing the
  number of labels (fuzzy terms).}

  \item{type.mf}{the type of the shape of the membership
  functions. See \code{\link{fuzzifier}}.}

  \item{type.tnorm}{the type of t-norm. See
  \code{\link{inference}}.}

  \item{type.snorm}{the type of s-norm. See
  \code{\link{inference}}.}

  \item{type.implication.func}{the type of implication
  function. See \code{\link{WM}}.}
}
\description{
  This is the internal function that implements the fuzzy
  rule-based classification system using Chi's technique
  (FRBCS.CHI). It is used to solve classification tasks.
  Users do not need to call it directly, but just use
  \code{\link{frbs.learn}} and \code{\link{predict}}. This
  method is suitable only for classification problems.
}
\details{
  This method was proposed by Zheru Chi, Hong Yan, and Tuan
  Pham that extends Wang and Mendel's method. The method
  consists of the following five steps: \itemize{ \item
  Step 1: Fuzzify the input space. \item Step 2: Generate
  fuzzy rules from given training data pairs. \item Step 3:
  Assign a degree to each rule. \item Step 4: Create a
  combined rule bank. \item Step 5: Determine the mapping
  by using a defuzzification method. }
}
\references{
  Z. Chi, H. Yan, T. Pham, "Fuzzy algorithms with
  applications to image processing and pattern
  recognition", World Scientific, Singapore (1996).
}
\seealso{
  \code{\link{FRBCS.eng}}, \code{\link{frbs.learn}}, and
  \code{\link{predict}}
}

