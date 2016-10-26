#' uhcsimXvalid
#'
#' \code{uhcsimXvalid} samples, randomly, locations from non-stratified
#' test data using k-fold cross validation.
#'
#' This function samples, randomly, locations from non-stratified test data and
#' returns an array of dimension nsims x nused_test x p (where p
#' is the number of predictors to be validated)
#'
#' @param data_sim The data set to be used in cross-validation
#' @param k_folds The number of bins to split data into
#' @param nsims The number of simulations (M) used to create the UHC plot.
#' @param nused The number of used locations in the data set.
#' @param navail The number of available locations in the data set (ie not used)
#' @param ncovariates The number of covariates or p in the model
#' @param model_form The model in the form: "y~x1+x2+...+xp" (including
#' quotations)
#' @param z_colnames The column names of the environmental characteristics to be
#' plotted with \code{uhcdensplot}. Takes the form: c("x1","x2",...,"xp")
#' @param xmat_colnames The column names of the environmental characteristics
#' that were used as predictors in the model_form. Takes the form:
#' c("x1","x2",...,"xp")
#'
#' @return An array of dimensions nsims x nused_test x p.
#'
#' @seealso Full archive of the data and code necessary to replicate the
#' manuscript at \url{http://doi.org/10.13020/D6T590}.
#'
#' @examples
#'
#' @export
uhcsimXvalid <- function(data_sim,
                         k_folds,
                         nsims,
                         nused,
                         navail,
                         ncovariates,
                         model_form,
                         z_colnames,
                         xmat_colnames){

  x_sim_choice <- array(NA,dim=c(nsims, nused, ncovariates))

  for(kk in 1:k_folds){
    # a) Fit model
    (glm_fit <- glm(model_form, family=binomial(), data=data_sim[data_sim$k!=kk,]))

    # Calculate the number of locations (used) in the "test" data bin
    #     This becomes the number of predicted "used" locations to draw later
    nused_test <- nrow(data_sim[data_sim$k==kk&data_sim$y==1,])

    # A vector or matrix of (used & available) environmental
    # characteristics in the test data set.
    z <- as.matrix(data_sim[data_sim$k==kk, z_colnames])
    # A matrix of predictor variables in the test data.
    xmat <- as.matrix(data_sim[data_sim$k==kk, xmat_colnames])

    # number of predictors in the model
    p <- length(coef(glm_fit))

    # b) Draw beta^i from mvn(hat{beta}, cov(\hat{\beta}))
    beta_hats <- MASS::mvrnorm(nsims,coef(glm_fit)[-1], vcov(glm_fit)[2:p,2:p])
    ntot_test <- nrow(xmat) # total number of test observations

    records <- 1:10
    records <- as.integer(records+(10*(kk-1)))
    for(mm in 1:nsims){
      # c) Calculate probability of selection based on $\beta^i$
      wx_pred <- exp(as.matrix(xmat)%*%beta_hats[mm,])
      # d) Predict which plots would be "used" from "test" data bin based on
      # probability of selection
      inds <- sample(1:ntot_test, nused_test, replace=FALSE, prob=wx_pred)
      # e) Record the covariates for those predicted used plots into the array
      x_sim_choice[mm,records,] <- z[inds,]
    }
  }

  return(x_sim_choice)
}
