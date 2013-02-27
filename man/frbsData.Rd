\docType{data}
\name{frbsData}
\alias{frbsData}
\title{Data set of the package}
\description{
  The package includes embedded versions of the
  Mackey-Glass chaotic time series and the Gas Furnance
  dataset.
}
\details{
  \bold{Mackey-Glass chaotic time series}

  The Mackey-Glass chaotic time series is defined by the
  following delayed differential equation:

  dx(t) / dt = (a * x(t - tau) / (1 + x(t - tau) ^ 10)) - b
  * x(t)

  For this dataset, we generated 1000 samples, with input
  parameters as follows: \itemize{ \item a = 0.2 \item b =
  0.1 \item tau = 17 \item x0 = 1.2 \item dt = 1 }

  The dataset is embedded in the following way:

  input variables: x(t - 18), x(t - 12), x(t - 6), x(t)

  output variable: x(t + 6)

  \bold{Gas Furnance dataset}

  The Gas Furnance dataset is taken from Box and Jenkins.
  It consists of 292 consecutive values of methane at time
  (t - 4), and the CO2 produced in a furnance at time (t -
  1) as input variables, with the produced CO2 at time (t)
  as an output variable. So, each training data point
  consists of [u(t - 4), y(t - 1), y(t)], where u is
  methane and y is CO2.
}
\references{
  Box, G. E. P., & Jenkins, G. M. "Time Series Analysis,
  forecasting and control", San Fransisco, CA: Holden Day
  (1970).

  Mackey, M., & Glass, L., "Oscillation and chaos in
  physiological control systems", Science, vol. 197, pp.
  287 - 289 (1977).
}
\keyword{data}
\keyword{datasets}

