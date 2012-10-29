\name{data.gen3d}
\alias{data.gen3d}
\title{A data generator}
\usage{
  data.gen3d(range.input, num.grid = 10)
}
\arguments{
  \item{range.input}{the range of the input variables, as a
  matrix(2 x n).}

  \item{num.grid}{a number representing the size of the
  grid on the plane.}
}
\value{
  the data
}
\description{
  The purpose of this function is to generate data, which
  contains two input variables and one output variable,
  automatically for all values on a plane.
}
\examples{
range.input <- matrix(c(0, 100, 0, 100), nrow=2)
num.grid <- 10
data.test <- data.gen3d(range.input, num.grid)
}

