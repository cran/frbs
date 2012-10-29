\name{frcs}
\alias{frcs}
\title{frcs model building}
\usage{
  frcs(range.data, data.train, label.inp, num.class,
    type.mf)
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

  \item{label.inp}{a matrix(1 x n) whose elements represent
  the number of labels (fuzzy terms), where n is the number
  of variables.}

  \item{num.class}{an integer number representing the
  number of labels (fuzzy terms).}

  \item{type.mf}{the type of the shape of the membership
  functions.}
}
\description{
  This is the internal function that implements the fuzzy
  rule-based classification system (frcs) model. Users do
  not need to call it directly, but just use
  \code{\link{frbs.learn}} and \code{\link{predict}}. This
  method is suitable only for classification problems.
}
\details{
  This method is adopted from Hisao Ishibuchi and Tomoharu
  Nakashima's paper and Oscar Cordon, Maria Jose del Jesus
  and Francisco Herrera's paper. The difference of a usual
  FRBS and this method is on the consequent part. The
  consequent part of frcs is a class with a certainty grade
  (i.e. rule weight). This method uses the techniques of
  Wang and Mendel for determining the antecedent part.
  Whereas on the consequent part, the value is associated
  to a class of the training data. After getting the fuzzy
  IF-THEN rules, we calculate a certainty grade (CF) for
  each rule. The prediction phase of this model is
  conducted by the internal function
  \code{\link{frcs.eng}}.
}
\references{
  Hisao Ishibuchi and Takashi Yamamoto, "Rule weight
  specification in fuzzy rule-based classification
  systems," IEEE Transactions on Fuzzy Systems, Vol. 13,
  No. 4 (2005).

  Hisao Ishibuchi and Tomoharu Nakashima, "Effect of rule
  weights in fuzzy rule-based classification systems", IEEE
  Transactions on Fuzzy Systems, Vol. 9, No. 4 (2001).

  Oscar Cordon, Maria Jose del Jesus and Francisco Herrera,
  "A proposal on reasoning methods in fuzzy rule-based
  classification systems", International Journal of
  Approximate Reasoning 20, pp. 21 - 45 (1999).

  Tomoharu N., Yasuyuki Y., Hisao Isibuchi, "Learning fuzzy
  IF-THEN rules for pattern classification with weighted
  training patterns," Proceedings of Joint 4th Conference
  of the European Society for Fuzzy Logic and Technology
  and the 11th Rencontres Francophones sur la Logique Floue
  et ses Applications, (2005).
}
\seealso{
  \code{\link{frcs.eng}}, \code{\link{frbs.learn}}, and
  \code{\link{predict}}
}

