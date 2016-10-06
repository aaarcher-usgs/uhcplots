#' uhcdiffdensplot
#'
#' \code{uhcdiffdensplot} plots an alternative UHC plot with a simulation
#' envelope for \eqn{f^U(z) - \hat{f}^u(z)}
#'
#' @inheritParams uhcdensplot
#'
#' @return A plot with the mean (black line) and upper and lower bounds of a
#' 95% simulation envelope (red dashed lines) of the difference
#' between the predicted density estimates and the density estimates for the
#' presence locations.
#'
#' @seealso Full archive of the data and code necessary to replicate the
#' manuscript at \url{http://doi.org/10.13020/D6T590}.
#'
#' @examples
#' # Simulate training data for the non-linear example
#' nonlinear.train <- uhcdatasimulator(nused = 100,
#'     navail = 10000,
#'     betas = c(2,-1),
#'     ntemp = 1000000,
#'     example = "non-linear")
#'
#' # Simulate test data for the non-linear example
#' nonlinear.test <- uhcdatasimulator(nused = 100,
#'     navail = 10000,
#'     betas = c(2,-1),
#'     ntemp = 1000000,
#'     example = "non-linear")
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
#'
#' # Get density estimates for quadratic model
#' denshats.correct <- uhcdenscalc(rand_sims = xhat.correct[,,1],
#'    dat = subset(nonlinear.test, y==1, select="temp"),
#'    avail = subset(nonlinear.test, y==0, select="temp"))
#'
#' # Get density estimates for linear (misspecified) model
#' denshats.misspec <- uhcdenscalc(rand_sims = xhat.misspec[,,1],
#'    dat = subset(nonlinear.test, y==1, select="temp"),
#'    avail = subset(nonlinear.test, y==0, select="temp"))
#'
#' # Create an alternate UHC plot for the quadratic model
#' uhcdiffdensplot(densdat = denshats.correct$densdat,
#'    densrand = denshats.correct$densrand)
#'
#' # Create an alternate UHC plot for the linear (misspecified) model
#' uhcdiffdensplot(densdat = denshats.misspec$densdat,
#'    densrand = denshats.misspec$densrand)
#' @export
uhcdiffdensplot <- function(densdat, densrand, xl=NULL){

  mean.f <- apply(densrand,2,mean)
  low.f <- apply(densrand,2,quantile, prob=0.025)
  up.f <- apply(densrand,2,quantile, prob=0.975)

  # Difference in density plot (f^u - hat(f)^u)
  yl2 <- range(c(densdat$y-low.f, densdat$y-up.f))
  if (is.null(xl)!= TRUE){
    plot(densdat$x, densdat$y-mean.f, ylim=yl2,
                              xlab="", type="l", ylab="", xlim=xl)
  }else{
      plot(densdat$x, densdat$y-mean.f, ylim=yl2, xlab="", type="l", ylab="")}

  lines(densdat$x, densdat$y-low.f, col="red", lty=2)
  lines(densdat$x, densdat$y-up.f, col="red", lty=2)
  abline(h=0)
}
