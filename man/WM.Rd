\name{WM}
\alias{WM}
\title{WM model building}
\usage{
  WM(data.train, num.labels, type.mf = "GAUSSIAN",
    type.tnorm = "PRODUCT",
    type.implication.func = "ZADEH",
    classification = FALSE, range.data = NULL)
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

  \item{type.mf}{the type of the membership function. See
  \code{\link{frbs.learn}}.}

  \item{type.tnorm}{a value which represents the type of
  t-norm. See \code{\link{inference}}.}

  \item{type.implication.func}{a value representing type of
  implication function. Let us consider a rule, a -> b,
  \itemize{ \item \code{DIENES_RESHER} means (b > 1 - a? b
  : 1 - a). \item \code{LUKASIEWICZ} means (b < a ? 1 - a +
  b : 1). \item \code{ZADEH} means (a < 0.5 || 1 - a > b ?
  1 - a : (a < b ? a : b)). \item \code{GOGUEN} means (a <
  b ? 1 : b / a). \item \code{GODEL} means (a <= b ? 1 :
  b). \item \code{SHARP} means (a <= b ? 1 : 0). \item
  \code{MIZUMOTO} means (1 - a + a * b). \item
  \code{DUBOIS_PRADE} means (b == 0 ? 1 - a : (a == 1 ? b :
  1)). \item \code{MIN} means (a < b ? a : b) }}

  \item{classification}{a boolean representing whether it
  is a classification problem or not.}

  \item{range.data}{a matrix representing interval of
  data.}
}
\description{
  This is the internal function that implements the model
  proposed by L.X. Wang and J.M. Mendel. It is used to
  solve regression task. Users do not need to call it
  directly, but just use \code{\link{frbs.learn}} and
  \code{\link{predict}}
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
  L. X. Wang and J.M. Mendel, "Generating fuzzy rule by
  learning from examples", IEEE Trans. Syst., Man, and
  Cybern., vol. 22, no. 6, pp. 1414 - 1427 (1992).
}
\seealso{
  \code{\link{frbs.learn}}, \code{\link{predict}} and
  \code{\link{frbs.eng}}.
}

