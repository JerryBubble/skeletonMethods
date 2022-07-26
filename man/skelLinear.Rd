% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skeletonRegression.R
\name{skelLinear}
\alias{skelLinear}
\title{Fitting the S-Lspline regressor}
\usage{
skelLinear(newnn, newpx, px, trainY, skeleton)
}
\arguments{
\item{newnn}{the matrix where each row records the indices of the two closest knots of the new test points.}

\item{newpx}{the projection proportion for the new test points.}

\item{px}{the projection proportion for the training points.}

\item{trainY}{the responses of the training points.}

\item{skeleton}{a list returned by the voronSkeleton method constructing the skeleton graph.}
}
\value{
A list with following components:
\describe{
\item{trainZ}{The transformed training data matrix.}
\item{testZ}{The transformed test data matrix.}
\item{model}{The lm object for the resulting regression model.}
\item{pred}{The vector of predicted responses for the test points.}
}
}
\description{
Fitting the S-Lspline regressor
}
