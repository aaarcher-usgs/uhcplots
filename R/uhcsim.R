#' uhcsim
#'
#' \code{uhcsim} samples, randomly, locations from non-stratified
#' test data.
#'
#' This function samples, randomly, locations from non-stratified test data and
#' returns an array of dimension nsims x nused_test x p (where p
#' is the number of predictors to be validated)
#'
#' @param nsims The number of simulations (M) used to create the UHC plot.
#' @param nused_test The number of used locations in the test data set.
#' @param xmat A matrix of predictor variables in the test data.
#' @param fit_rsf The fitted logistic regression model object
#' @param z A vector or matrix of (used & available) environmental
#' characteristics in the test data set.
#'
#' @return An array of dimensions nsims x nused_test x p.
#'
#' @seealso Full archive of the data and code necessary to replicate the
#' manuscript at \url{http://doi.org/10.13020/D6T590}.
#'
#' @examples
#' # Simulate training data for the non-linear example
#' nonlinear.train <- uhcdatasimulator(nused = 100,
#'    navail = 10000,
#'    betas = c(2,-1),
#'    ntemp = 1000000,
#'    example = "non-linear")
#'
#' # Simulate test data for the non-linear example
#' nonlinear.test <- uhcdatasimulator(nused = 100,
#'    navail = 10000,
#'    betas = c(2,-1),
#'    ntemp = 1000000,
#'    example = "non-linear")
#'
#' # Fit GLM with quadratic relationship
#' train.correct <- glm(y~temp + I(temp^2),
#'    family = binomial,
#'    data = nonlinear.train)
#'
#' # Fit GLM with linear (misspecified) relationship
#' train.misspec <- glm(y~temp,
#'    family = binomial,
#'    data = nonlinear.train)
#'
#' # Simulate data for quadratic model
#' xhat.correct <- uhcsim(nsims = 1000,
#'    nused_test = 100,
#'    xmat = model.matrix(y~temp + I(temp^2), data = nonlinear.test)[,-1],
#'    fit_rsf = train.correct,
#'    z = as.matrix(nonlinear.test[,"temp"]))
#'
#' # Simulate data for linear (misspecified) model
#' xhat.misspec <- uhcsim(nsims = 1000,
#'    nused_test = 100,
#'    xmat = as.matrix(model.matrix(y~temp, data = nonlinear.test)[,2]),
#'    fit_rsf = train.misspec,
#'    z = as.matrix(nonlinear.test[,"temp"]))
uhcsim <- function(nsims, nused_test, xmat, fit_rsf, z){
  # array to store chosen x's for each simulation
  x_sim_choice <- array(NA,dim=c(nsims,nused_test,ncol(z)))
  p <- length(coef(fit_rsf)) # number of predictors in the model

  # new beta^ for each simulation
  beta.hats <- MASS::mvrnorm(nsims,coef(fit_rsf)[-1], vcov(fit_rsf)[2:p,2:p])
  ntot.test <- nrow(xmat) # total number of test observations
  z <- as.matrix(z) # z has to be a matrix for the code to work

  for(i in 1:nsims){
    # probability of choosing each location
    wx.pred <- exp(as.matrix(xmat)%*%beta.hats[i,])
    # choices
    inds <- sample(1:ntot.test, nused_test, replace=FALSE, prob=wx.pred)
    # covariates at chosen locations
    x_sim_choice[i,,] <- z[inds,]
  }
  return(x_sim_choice)
}
