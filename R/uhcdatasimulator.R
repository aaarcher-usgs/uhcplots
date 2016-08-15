#' uhcdatasimulator
#'
#' \code{uhcdatasimulator} simulates the data used in Fieberg et al. (In Review)
#'
#' This is a function that creates a dataframe based on the example and scenario
#' chosen from the manuscript (Fieberg et al. In Review). In the first example,
#' ("missing predictor") the distribution of a species is related to elevation
#' and precipitation. In the second example ("non-linear"), the distribution of
#' a species is non-linearly related to temperature.
#'
#' @param nused The number of used locations in training/test data set
#' @param navail The number of background locations in training/test data set
#' @param betas The vector of length 2 for *true* probability of use
#' @param corx The correlation between elevation and precipitation in
#' training/test dataset. For missing predictor example only.
#' @param ntemp A large number of available points.
#' @param scenario The specific scenario as described in the manuscript.
#' Options include"corr00", "corrNN", "corrPN", and "T".
#' @param example The name of the example. Options include "missing predictor"
#' or "non-linear".
#'
#' @return A dataframe of simulated data
#'
#' @examples
#' # Simulate data for the non-linear example
#' new.data <- uhcdatasimulator(nused = 100,
#'    navail = 10000,
#'    betas = c(2,-1),
#'    ntemp = 1000000,
#'    scenario = "T",
#'    example = "non-linear")
uhcdatasimulator <- function(nused, navail,
                             betas, corx, ntemp,
                             scenario, example){
  ### Missing predictor example
  if (example=="missing predictor"){
    # ntot = total number of locations in data set
    ntot <- nused + navail

    # Create elevation and precipitation data
    covx <- corx*4
    xdat <- MASS::mvrnorm(navail, mu=c(0,0),
                          Sigma=matrix(c(4,  covx,  covx, 4),
                                       ncol=2, byrow=TRUE))
    elev <- xdat[,1] # elevation data
    precip <- xdat[,2] # precipitation data

    # Generate *used* locations by selecting these points from from a large #
    #    (ntemp) of available points, and assemble data
    tempdat1 <- MASS::mvrnorm(ntemp, mu=c(0,0),
                        Sigma=matrix(c(4,covx, covx, 4),
                                     ncol=2, byrow=TRUE))
    elevtemp <- tempdat1[,1] # available elev
    preciptemp <- tempdat1[,2] # available precip

    # used in data based on *true* probability (wx.true)
    wx.true <- exp(elevtemp*betas[1]+preciptemp*betas[2])
    inds <- sample(1:ntemp, size=nused, replace=TRUE,
                         prob=wx.true)

    # Compile data for scenario
    sim.data <- data.frame(y=c(rep(0, navail), rep(1, nused)),
                                  elev=c(elev, elevtemp[inds]),
                                  precip=c(precip, preciptemp[inds]))
  }
  ### Non-linear example
  if (example=="non-linear"){
    # Generate *available* covariates in training and test data sets
    #       (assume these distributions are the same in the test and
    #         training data sets)
    temp <- rnorm(navail, 0, 2)

    # Generate *used* locations by selecting from a large number
    #     (ntemp, below) of available points, and assemble data
    temptemp <- rnorm(ntemp, 0, 2) # available points to choose from

    # used in data based on *true* probability (wx.true)
    wx.true <- exp(temptemp*betas[1]+temptemp^2*betas[2])
    inds <- sample(1:ntemp, size=nused,
                         replace=TRUE, prob=wx.true)


    # Compile data
    sim.data <- data.frame(y=c(rep(0, navail), rep(1, nused)),
                             temp=c(temp, temptemp[inds]))

  }
  # Return files
  return(sim.data)
}
#' spun with ezspin("uhcplots/functions/uhcdatasimulator.R", out_dir = "uhcplots/output", fig_dir="uhcplots/figures", keep_md=F)
