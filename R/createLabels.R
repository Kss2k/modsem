

createLabelCov <- function(x, y) {
  paste("Cov", x, y, sep = "_")
}



createLabelVar <- function(x) {
  paste("Var", x, sep = "_")
}



createLabelLambda <- function(ind, latent) {
  paste("lambda", ind, latent, sep = "_")
}


createLabelLambdaSquared <- function(ind, latent) {
  paste0(createLabelLambda(ind, latent), " ^ 2")
}
