\name{WM}
\alias{WM}
\title{WM model building}
\usage{
  WM(range.data, data.train, num.labels = 5, type.mf = 3,
    classification = FALSE)
}
\arguments{
  \item{range.data}{a matrix(2 x n) containing the range of
  the data, where n is the number of variables, and first
  and second rows are the minimum and maximum values,
  respectively.}

  \item{data.train}{a matrix(m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables; the last column is the
  output variable.}

  \item{num.labels}{a matrix(1 x n), whose elements
  represent the number of labels (fuzzy terms); n is the
  number of variables.}

  \item{type.mf}{the type of the membership function. See
  \code{\link{frbs.learn}}.}

  \item{classification}{a boolean representing whether it
  is a classification problem or not.}
}
\description{
  This is the internal function that implements the model
  proposed by L.X. Wang and J.M. Mendel. Users do not need
  to call it directly, but just use
  \code{\link{frbs.learn}} and \code{\link{predict}}
}
\details{
  The fuzzy rule-based system for learning from L.X. Wang
  and J.M. Mendel's paper is implemented in this function.
  For the learning process, there are three stages: Step 1
  divides the input and output spaces of the given
  numerical data into fuzzy regions. Step 2 generates fuzzy
  IF-THEN rules from the training data. Step 3 determines a
  degree of each rule. In the prediction phase, there are
  four steps: fuzzification, checking the rules, inference,
  and defuzzification.
}
\references{
  L. X. Wang, "Fuzzy systems are universal approximators,"
  in Proc. IEEE Int. Conf. Fuzzy Systems, San Diego, CA
  (1992).

  L. X. Wang and J.M. Mendel, "Generating fuzzy rule by
  learning from examples", in Proc. 6th Int. Symp.
  Intelligent Control (Washington, DC), 1991, pp. 263-268;
  also IEEE Trans. Syst., Man, Cybern., vol. 22, No. 6,
  (1992).

  L. X. Wang and J. M. Mendel, "Fuzzy basis function,
  universal approximation, and orthogonal least squares
  learning," IEEE Int. Conf. Neural Network, vol. 3 no. 5,
  pp. 807 - 814 (1992).
}
\seealso{
  \code{\link{frbs.learn}}, \code{\link{predict}} and
  \code{\link{frbs.eng}}.
}

