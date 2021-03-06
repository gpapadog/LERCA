% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenYgivenXC_function.R
\name{GenYgivenXC}
\alias{GenYgivenXC}
\title{Generating the outcome given exposure and covariates.}
\usage{
GenYgivenXC(dataset, out_coef, bY, bYX, Ysd, XY_function, XY_spec = NULL)
}
\arguments{
\item{dataset}{A data frame containing X, C1, C2, ..., and E.}

\item{out_coef}{A matrix where the ij element is the coefficient of
covariate Ci in the outcome model of experiment j.}

\item{bY}{Vector of length equal to the number of experiments including the
intercept of the outcome model in each experiment.}

\item{bYX}{If XY_function is set to linear, bYX includes the coefficient of
the exposure within each experiment and is of length num_exper. If
XY_function is set to other, bYX is of length 1 and corresponds to the
coefficient in front of the exposure term.}

\item{Ysd}{Vector of length number of experiments. Standard deviation of
outcome model residual.}

\item{XY_function}{String specifying whether the XY relationship is piece-
wise linear (set 'linear'), or a continuous function supplied by the XY_spec
arguement (set 'other').}

\item{XY_spec}{Needs to be specified if XY_function is set to 'other'. It is
the function that specifies the true ER relationship. Defaults to NULL.}
}
\description{
In the simulated data, we specify the outcome model and residual variance
for every experiment and the outcome is generated.
}
