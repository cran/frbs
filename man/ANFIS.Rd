\name{ANFIS}
\alias{ANFIS}
\title{ANFIS model building}
\usage{
  ANFIS(data.train, num.labels, max.iter = 10,
    step.size = 0.01, type.tnorm = "MIN",
    type.snorm = "MAX", type.implication.func = "ZADEH")
}
\arguments{
  \item{data.train}{a matrix(m x n) of normalized data for
  the training process, where m is the number of instances
  and n is the number of variables; the last column is the
  output variable.  Note the data must be normalized
  between 0 and 1.}

  \item{num.labels}{a matrix(1 x n), whose elements
  represent the number of labels (fuzzy terms); n is the
  number of variables.}

  \item{max.iter}{the maximal number of iterations.}

  \item{step.size}{a real number between 0 and 1
  representing the step size of the gradient descent.}

  \item{type.tnorm}{the type of t-norm. For more detail,
  please have a look at \code{\link{inference}}.}

  \item{type.snorm}{the type of s-norm. For more detail,
  please have a look at \code{\link{inference}}.}

  \item{type.implication.func}{a value representing type of
  implication function.  For more detail, please have a
  look at \code{\link{WM}}.}
}
\description{
  This is the internal function that implements the
  adaptive-network-based fuzzy inference system (ANFIS). It
  is used to solve regression tasks. Users do not need to
  call it directly, but just use \code{\link{frbs.learn}}
  and \code{\link{predict}}.
}
\details{
  This method was proposed by Jyh-Shing and Roger Jang. It
  uses the Takagi Sugeno Kang model on the consequent part
  of the fuzzy IF-THEN rules. The ANFIS architecture
  consists of two processes, the forward and the backward
  stage. The forward stage has five layers as follows:
  \itemize{ \item Layer 1: The fuzzification process which
  transforms crisp values into fuzzy terms using the
  Gaussian function as the shape of the membership
  function. \item Layer 2: The inference stage using the
  t-norm operator (the AND operator). \item Layer 3:
  Calculating the ratio of the strengths of the rules.
  \item Layer 4: Calculating the consequent parameters.
  \item Layer 5: Calculating the overall output as the sum
  of all incoming signals. } The backward stage is a
  process of parameter learning. In this step, the least
  squares method is used in order to obtain the parameters,
  which are coefficients of linear equations on the
  consequent part, and mean and variance on the antecedent
  part.
}
\references{
  Jyh-Shing and Roger Jang, "ANFIS: adaptive-network-based
  fuzzy inference system", IEEE Transactions on Systems,
  Man, and Cybernetics, vol. 23, no. 3, pp. 665 - 685
  (1993).

  Jyh-Shing Roger Jang, et.al., "Neuro-fuzzy and soft
  computing: a computational approach to learning and
  machine intelligence", Prentice-Hall, Inc (1997).
}
\seealso{
  \code{\link{ANFIS.update}}, \code{\link{frbs.learn}}, and
  \code{\link{predict}}
}

