#' uhcdensplot
#'
#' \code{uhcdensplot} creates UHC plots
#'
#' This function plots the density of the environmental
#' characteristics at the observed locations in the test data set, $f^u(z)$,
#' along with a simulation envelope for $f^U(z)$ created by randomly choosing
#' locations in the test data using the \code{uhcsim} or \code{uhcsimstrat}
#' function.
#'
#' @param densdat The kernel density estimates of observed points in the test
#' data set
#' @param densrand The predicted kernel density estimates of simulated data
#' points
#' @param includeAvail An indicator determining whether the distribution of
#' the available locations should be drawn on the plot
#' @param densavail The kernel density estimates of available points in the
#' test data set
#' @param xl The x-axis limits (can be user supplied)
#' @param yl The y-axis limits (can be user supplied)
#' @param includeLegend An indicator determining whether the legend should
#' be displayed.
#'
#' @return A plot with a 95% simulation envelope (grey band) of the predicted
#' density estimates for simulated data points and kernel density estimates of
#' the environmental characteristics associated with the observed locations
#' (red dashed line) and the available locations (black line).
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
#' poop
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
#' # Create a UHC plot for the quadratic model
#' uhcdensplot(densdat = denshats.correct$densdat,
#'    densrand = denshats.correct$densrand,
#'    includeAvail = TRUE,
#'    densavail = denshats.correct$densavail,
#'    includeLegend = TRUE)
#'
#' # Create a UHC plot for the linear (misspecified) model
#' uhcdensplot(densdat = denshats.misspec$densdat,
#'    densrand = denshats.misspec$densrand,
#'    includeAvail = TRUE,
#'    densavail = denshats.misspec$densavail,
#'    includeLegend = TRUE)
uhcdensplot <- function(densdat, densrand, includeAvail=F, densavail=NULL, xl=NULL,
                        yl=NULL, includeLegend=T){
  # combine to get reasonable axis limits
  alldens <- c(densrand, densdat$y, densavail$y)
  if (is.null(yl)){
  yl <- c(min(alldens), max(alldens))
  }

  # mean predicted density
  mean.f <- apply(densrand,2,mean, na.rm=TRUE)
  #lower sim envelope
  low.f <- apply(densrand,2,quantile, prob=0.025, na.rm=TRUE)
  # upper sim envelope
  up.f <- apply(densrand,2,quantile, prob=0.975, nar.rm=TRUE)
  if (is.null(xl)!= TRUE){
    plot(densdat$x, densdat$y, ylim=yl, xlab="",
                               type="n", ylab="Density", xlim=xl)
  }else{
      plot(densdat$x, densdat$y, ylim=yl, xlab="", type="n", ylab="")
    }

  # Observed and predicted density of f^u
  polygon(x = c(densdat$x, rev(densdat$x)),
          y = c(up.f, rev(low.f)), col="gray", border="gray")
  lines(densdat$x, densdat$y, ylim=yl, xlab="",
        type="l", ylab="", col="red", lty=2, lwd=2)

  # Include available if includeAvail=TRUE
  if (includeAvail==TRUE){
    lines(densavail$x, densavail$y, lwd=2)
  }

  if (includeAvail==TRUE & includeLegend==TRUE){
    legend(-5.2, max(alldens),
           c("Available", "Used", "Predicted"),
           lty = c(1,2,1),
           col = c("black", "red", "grey"),
           bty = "n",
           lwd = c(1,1,5))
  }

  if (includeAvail==F & includeLegend==TRUE){
    legend(-5.2, max(alldens),
           c("Used", "Predicted"),
           lty = c(2,1),
           col = c("red", "grey"),
           bty = "n",
           lwd = c(1,5))
  }
}
