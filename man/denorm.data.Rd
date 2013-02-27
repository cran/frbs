\name{denorm.data}
\alias{denorm.data}
\title{The data de-normalization}
\usage{
  denorm.data(dt.norm, range.data, min.scale = 0,
    max.scale = 1)
}
\arguments{
  \item{dt.norm}{a matrix(n x m) of the normalized data.}

  \item{range.data}{a matrix(2 x n) containing the range of
  the data, where n is the number of variables, and first
  and second rows are the minimum and maximum value,
  respectively.}

  \item{min.scale}{the minimum value within normalization.}

  \item{max.scale}{the maximum value within normalization.}
}
\value{
  the real-valued data
}
\description{
  This function is to transform from normalized data into
  real-valued data.
}
\seealso{
  \code{\link{norm.data}}
}

