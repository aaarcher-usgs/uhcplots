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
#' @param model_form The model in the form: "y~x1+x2+...+xp" (including quotations)
#' @param z_colnames The column names of the environmental characteristics to
#' be plotted with \code{uhcdensplot}. Takes the form: c("x1","x2",...,"xp")
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
#' # Simulate data using uhcplots::uhcdatasimulator
#' data_sim <- uhcdatasimulator(nused = 100,
#'                             navail = 10000,
#'                             betas = c(0.5, -1) ,
#'                             corx = 0,
#'                             ntemp = 100000,
#'                             example = "missing predictor")
#'
#' #Set up parameters
#'
#' # Number of used and available locations in dataset
#' used <- data_sim[data_sim$y==1,];  nused <- nrow(used)
#' avail <- data_sim[data_sim$y==0,]; navail <- nrow(avail)
#'
#' ## k-fold cross validation
#'
#' #**Step 1** Partition data into k-fold bins
#' k_folds <- 10
#' nused_k <- nused/k_folds
#' navail_k <- navail/k_folds
#' data_sim$k <- dismo::kfold(data_sim, k = k_folds, by = data_sim$y)
#' table(data_sim$k)
#'
#' #**Step 2** Fit a model (in this case, missing precipitation)
#'
#' # Set up some parameters
#' nsims <- 1000 # M from paper;
#' ncovariates <- 2;
#' model_form_missingP <- "y~elev"
#' z_colnames <- c("elev","precip")
#' xmat_colnames_missingP <- c("elev")
#'
#' # Create an array to store chosen location covariates for each simulation
#' # and each k-fold "test" data bin
#' x_sim_choice_missingP <- array(NA, dim = c(nsims, nused, ncovariates))

#' # Get predicted distribution of covariates for each fold, k:
#' x_sims_choice_missingP <- uhcplots::uhcsimXvalid(data_sim = data_sim,
#'                                                 k_folds = k_folds,
#'                                                 nsims = nsims,
#'                                                 nused = nused,
#'                                                 navail = navail,
#'                                                 ncovariates = ncovariates,
#'                                                 model_form = model_form_missingP,
#'                                                 z_colnames = z_colnames,
#'                                                 xmat_colnames = xmat_colnames_missingP)
#' # **Step 3** Density calculations
#' denscalc_elev_missingP <- uhcdenscalc(rand_sims = x_sims_choice_missingP[,,1],
#'                                      dat = used[,"elev"],
#'                                      avail = avail[,"elev"],
#'                                      gridsize = 500)
#' denscalc_prec_missingP <- uhcdenscalc(rand_sims = x_sims_choice_missingP[,,2],
#'                                      dat = used[,"precip"],
#'                                      avail = avail[,"precip"],
#'                                      gridsize = 500)
#' # **Step 4** Create UHC plots for each covariate
#' opar<-par()
#' layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4,1))
#' par(mai=rep(0.5, 4))
#' uhcdensplot(densdat = denscalc_elev_missingP$densdat,
#'            densrand = denscalc_elev_missingP$densrand,
#'            includeAvail = T,
#'            densavail = denscalc_elev_missingP$densavail,includeLegend = F)
#' mtext(outer=F, side=1, line = 3, "Elevation")
#'
#' uhcdensplot(densdat = denscalc_prec_missingP$densdat,
#'            densrand = denscalc_prec_missingP$densrand,
#'            includeAvail = T,
#'            densavail = denscalc_prec_missingP$densavail,
#'            includeLegend = F)
#' mtext(outer=F, side=1, line= 3, "Precipitation")
#' par(mai=c(0,0,0,0))
#' plot.new()
#' legend("center",c("Available", "Used", "Predicted"), lty = c(1,2,1), col = c("black", "red", "grey"), bty = "n", lwd = c(1,1,5))
#' par<-opar
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

    records <- 1:nused_test
    records <- as.integer(records+(nused_test*(kk-1)))
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
