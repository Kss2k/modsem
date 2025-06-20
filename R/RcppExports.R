# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

calcSESimpleSlopes <- function(X, V) {
    .Call(`_modsem_calcSESimpleSlopes`, X, V)
}

muLmsCpp <- function(model, z) {
    .Call(`_modsem_muLmsCpp`, model, z)
}

sigmaLmsCpp <- function(model, z) {
    .Call(`_modsem_sigmaLmsCpp`, model, z)
}

completeLogLikLmsCpp <- function(modelR, P, quad) {
    .Call(`_modsem_completeLogLikLmsCpp`, modelR, P, quad)
}

gradLogLikLmsCpp <- function(modelR, P, block, row, col, eps = 1e-6) {
    .Call(`_modsem_gradLogLikLmsCpp`, modelR, P, block, row, col, eps)
}

muQmlCpp <- function(m, t) {
    .Call(`_modsem_muQmlCpp`, m, t)
}

sigmaQmlCpp <- function(m, t) {
    .Call(`_modsem_sigmaQmlCpp`, m, t)
}

calcKronXi <- function(m, t) {
    .Call(`_modsem_calcKronXi`, m, t)
}

calcBinvCpp <- function(m, t) {
    .Call(`_modsem_calcBinvCpp`, m, t)
}

dnormCpp <- function(x, mu, sigma) {
    .Call(`_modsem_dnormCpp`, x, mu, sigma)
}

varZCpp <- function(Omega, Sigma1, numEta) {
    .Call(`_modsem_varZCpp`, Omega, Sigma1, numEta)
}

#' Multiply indicators 
#' @param df A data DataFrame
#' @return A NumericVector
#' @export
multiplyIndicatorsCpp <- function(df) {
    .Call(`_modsem_multiplyIndicatorsCpp`, df)
}

rep_dmvnorm <- function(x, expected, sigma, t) {
    .Call(`_modsem_rep_dmvnorm`, x, expected, sigma, t)
}

dmvnrm_arma_mc <- function(x, mean, sigma, logd = TRUE) {
    .Call(`_modsem_dmvnrm_arma_mc`, x, mean, sigma, logd)
}

totalDmvnWeightedCpp <- function(mu, sigma, nu, S, tgamma, n, d) {
    .Call(`_modsem_totalDmvnWeightedCpp`, mu, sigma, nu, S, tgamma, n, d)
}

tracePathsNumericCpp <- function(x, y, parTable, maxlen = 100L) {
    .Call(`_modsem_tracePathsNumericCpp`, x, y, parTable, maxlen)
}

tracePathsCharacterCpp <- function(x, y, parTable, paramCol = "mod", maxlen = 100L) {
    .Call(`_modsem_tracePathsCharacterCpp`, x, y, parTable, paramCol, maxlen)
}

