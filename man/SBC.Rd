\name{SBC}
\alias{SBC}
\title{The subtractive clustering and fuzzy c-means (SBC) model building}
\usage{
  SBC(data.train, range.data.ori, r.a = 0.5,
    eps.high = 0.5, eps.low = 0.15)
}
\arguments{
  \item{data.train}{a matrix(m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables; the last column is the
  output variable.}

  \item{range.data.ori}{a matrix(2 x n) containing the
  range of the data, where n is the number of variables,
  and first and second rows are the minimum and maximum
  value, respectively.}

  \item{r.a}{the radius defining a neighborhood.}

  \item{eps.high}{an upper threshold value.}

  \item{eps.low}{a lower threshold value.}
}
\description{
  This is the internal function that implements a
  combination of the subtractive clustering method and
  fuzzy c-means. Users do not need to call it directly, but
  just use \code{\link{frbs.learn}} and
  \code{\link{predict}}
}
\details{
  This method was proposed by Stephen Chiu and J.C. Bezdek.
  For generating the rules in the learning phase, the
  subtractive clustering method is used to obtain the
  cluster centers. Subtractive clustering (SBC) is an
  extension of Yager and Filev's mountain method. SBC
  considers each data point as a potential cluster center
  by determining the potential of a data point as a
  function of its distances to all the other data points. A
  data point has a high potential value if that data point
  has many nearby neighbors. The highest potential is
  chosen as the cluster center and then the potential of
  each data point will be updated. The process of
  determining new clusters and updating potentials repeats
  until the remaining potential of all data points falls
  below some fraction of the potential of the first cluster
  center. After getting all the cluster centers from
  subtractive clustering, the cluster centers are optimized
  by fuzzy c-means.
}
\references{
  Bezdek, J.C., "Pattern recognition with fuzzy objective
  function algorithms", Plenum Press, New York (1981).

  Nikhil R. Pal, James C. Bezdek, and Richard J. Hathaway,
  "Sequential competitive learning and the fuzzy c-means
  clustering algorithms," Neural Networks, Vol. 9, No. 5,
  pp. 787-796 (1996).

  R. Yager and D. Filev, "Generation of fuzzy rules by
  mountain clustering," J. of Intelligent and Fuzzy
  Systems, vol. 2, no. 3, pp. 209 - 219 (1994).

  Stephen Chiu, "Method and software for extracting fuzzy
  classification rules by subtractive clustering", Fuzzy
  Information Processing Society, NAFIPS, pp. 461 - 465
  (1996).
}
\seealso{
  \code{\link{SBC.test}}, \code{\link{frbs.learn}}, and
  \code{\link{predict}}
}

