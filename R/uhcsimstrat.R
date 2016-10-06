#' uhcsimstrat
#'
#' \code{uhcsimstrat} samples, randomly, locations from stratified
#' test data.
#'
#' This function samples, randomly, locations from stratified
#' test data and returns an array of dimension nsims x nused_test x p (where p
#' is the number of predictors to be validated)
#'
#' @param stratum The stratum identifiers associated with each point in
#' xmat and z
#' @param fit_ssf The fitted step-selection function object
#' @inheritParams uhcsim
#'
#' @return An array of dimensions nsims x nused_test x p.
#'
#' @seealso Full archive of the data and code necessary to replicate the
#' manuscript at \url{http://doi.org/10.13020/D6T590}.
#'
#' @examples
#' # Load example moose data
#' mdat <- moose12687
#'
#' # Split into training and test datasets
#' mdat.train <- subset(mdat, year==2013)
#' mdat.test <- subset(mdat, year==2014)
#'
#' # Fit step-selection function with all covariates (Full model)
#' ssf.train.full <- survival::clogit(presence ~ decid50 + mixed50 +
#'    conif50 + treedwet50 + step + strata(stratum),
#'    data = mdat.train)
#'
#' # Fit step-selection function with only mixed50 (Reduced model)
#' ssf.train.reduc <- survival::clogit(presence ~ mixed50 +
#'    step + strata(stratum),
#'    data = mdat.train)
#'
#' # Create design matrix from test data for Full model SSF
#' design.mat.full <- model.matrix(~decid50 + mixed50 +
#'   conif50 + treedwet50 + step -1,
#'   data = mdat.test)
#'
#' # Create design matrix from test data for Reduced model SSF
#' design.mat.reduc <- model.matrix(~mixed50 + step -1, data = mdat.test)
#'
#' # Create design matrix for covariates for z (matrix of used & available
#' # environmental characteristics in the test data set.)
#' z <- model.matrix(~decid50 + mixed50 + conif50 + treedwet50 + step -1,
#'    data = mdat.test)[,-5]
#'
#' # Simulate data for Full model SSF
#' xhat.full <- uhcsimstrat(nsims = 1000,
#'   xmat = design.mat.full,
#'   stratum = mdat.test$stratum,
#'   fit_ssf = ssf.train.full,
#'   z = z)
#'
#' # Simulate data for Reduced model SSF
#' xhat.reduc <- uhcsimstrat(nsims = 1000,
#'   xmat = design.mat.reduc,
#'   stratum = mdat.test$stratum,
#'   fit_ssf = ssf.train.reduc,
#'   z = z)
#' @export
uhcsimstrat <- function(nsims, xmat, stratum, fit_ssf, z){

  ustrat <- unique(stratum) # unique strata
  nstrat <- length(ustrat) # number of strata

  # array to store chosen x's for each simulation
  x_sim_choice <- array(NA, dim = c(nsims, nstrat, ncol(z)))
  # new beta^ for each simulation
  beta.hats <- MASS::mvrnorm(nsims, coef(fit_ssf), vcov(fit_ssf))
  ntot.test <- nrow(xmat) # total number of test observations

  dimxmat <- ncol(xmat)
  tempdat <- as.data.frame(cbind(1:nrow(xmat),stratum, z))
  names(tempdat)[1] <- "inds"

  for(i in 1:nsims){
    # Use new beta^ vector to generate predicted values for test data
    lp.hat.s <- as.matrix(xmat)%*%beta.hats[i,]
    tempdat$wx.s <- exp(lp.hat.s)
    tempdat.g <- tempdat%>%group_by(stratum)
    x_sim_choice[i,,] <- z[sample_n(tempdat.g,1, weight=wx.s)$inds,] # choices
  }
  return(x_sim_choice)
}

