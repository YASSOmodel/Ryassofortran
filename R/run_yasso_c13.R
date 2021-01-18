#' Run the YASSO model
#'
#' \code{run_yasso()} runs the YASSO model and returns simulated soil carbon.
#'
#' \code{run_yasso()} wraps the Fortran90-release of the soil carbon model
#' YASSO15 into a simple R-function. The function is a convenient way to call
#' the Fortran-release.
#'
#' The function provides YASSO with the initial soil carbon values in the vector
#' \code{init} and runs the model one time step at a time. The simulated
#' carbon of the current time step is used as the initial value of the next time
#' step. The model runs until it has looped over all the time steps.
#'
#' @param par A numeric vector of YASSO parameters.
#' @param n_runs Input data. Refer to \code{\link{sample_data_run}} for now.
#' @param time -||-
#' @param temp -||-
#' @param prec -||-
#' @param init -||-
#' @param litter -||-
#' @param wsize -||-
#' @param leac -||-
#' @param sspred Optional integer, should steady state mode be used (1 = yes).
#'
#' @return A matrix containing the initial soil carbon on the first row and
#'   simulated soil carbon on the following rows.
#' @export
#'
#' @examples
#' soil_c <- run_yasso(
#'  par = sample_parameters,
#'  n_runs = sample_data_run$n_runs,
#'  time = sample_data_run$time,
#'  temp = sample_data_run$temp,
#'  prec = sample_data_run$prec,
#'  init = sample_data_run$init,
#'  litter = sample_data_run$litter,
#'  wsize = sample_data_run$wsize,
#'  leac = sample_data_run$leac
#' )

run_yasso_c13 <- function(par, par_c13, n_runs, time, temp, prec, soil_c, init, litter, wsize, leac,
	      	 	  sspred = 0L) {

  # Typeset parameters
  par <- as.double(par)
  par_c13 <- as.double(par_c13)
  
  # Initialize, typeset an array for the results
  soil_c13 <- matrix(rep(0, len = (n_runs) * 5), nrow = n_runs)
  soil_c13 <- rbind(init, soil_c13)
  soil_c13 <- unname(as.matrix(soil_c13))

  # Call the fortran model
  xx <- .Fortran(
    "runyasso_c13",
    par = par,
    par_c13 = par_c13,
    n_runs = n_runs,
    time = time,
    temp = temp,
    prec = prec,
    litter = litter,
    wsize = wsize,
    leac = leac,
    soil_c = soil_c,
    soil_c13 = soil_c13,
    sspred = sspred
  )

  # Return simulated soil carbon
  return(xx$soil_c13)
}