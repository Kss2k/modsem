% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{standardized_estimates}
\alias{standardized_estimates}
\title{Get standardized estimates}
\usage{
standardized_estimates(object, ...)
}
\arguments{
\item{object}{An object of class `modsem_da`,  `modsem_mplus`, 
or a parTable of class `data.frame`}

\item{...}{Additional arguments passed to other functions}
}
\description{
Get standardized estimates
}
\details{
for `modsem_lms`, `modsem_qml` and `modsem_mplus` objects, 
the interaction term is not standardized such that var(xz) = 1. 
The interaction term is not an actual variable in the model, meaning that it does not 
have a variance. It must therefore be calculated from the other parameters in the model.
Assuming normality and zero-means the variance is calculated as 
`var(xz) = var(x) * var(z) + cov(x, z)^2`. Thus setting the variance of the interaction 
term to 1, would only be 'correct' if the correlation between x and z is zero.
This means that the standardized estimates for the interaction term will 
be different from those using lavaan, since there the interaction term is an 
actual latent variable in the model, with a standardized variance of 1.
}
