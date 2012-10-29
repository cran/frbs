\name{norm.data}
\alias{norm.data}
\title{The data normalization}
\usage{
  norm.data(dt.ori, range.data, min.scale = -1,
    max.scale = 1)
}
\arguments{
  \item{dt.ori}{a matrix(n x m) of the original data.}

  \item{range.data}{a matrix(2 x n) containing the range of
  the data, where n is the number of variables, and first
  and second rows are the minimum and maximum value,
  respectively.}

  \item{min.scale}{the minimum value within normalization.}

  \item{max.scale}{the maximum value within normalization.}
}
\description{
  This function is to transform from real-valued data into
  normalized data.
}
\seealso{
  \code{\link{denorm.data}}
}

