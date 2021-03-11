#' Test data generator
#'
#' Generate test data in order to test the filtering functions.
#' To a signal function is added random noise contribution.
#' V0.1 = noise is assumed gaussian
#'
#' @param amplitude amplitude of the signal, default = 1
#' @param f frequency of the sinusoidal signal, default = 100
#' @param npoints number of points of the time serie
#' @param type type of signal, default = sinusoidal. Available types: sinusoidal, gaussian
#' @param x0 signal position for gaussian type. Default = 0
#' @param noise_contribution percentage pointing the maximum wanted signal/noise ratio. Default = 10
#'
#' @return A time serie with added random noise.
#'
#' @examples
#' test_data()
#'
#' @importFrom stats runif
#'
#' @export

test_data <- function(amplitude = 1, f = 100, npoints = 1000, type = "sinusoidal", x0 = 0, noise_contribution = 100){
  if (type == "sinusoidal") {
    t = seq(0, (1/f)*10, 1/((npoints/10)*f))
    time_serie = amplitude*sin(2*pi*f*t)
    noise = runif(npoints, min(time_serie)*(noise_contribution)/100, max(time_serie)*(noise_contribution)/100)
    return(time_serie[1:(length(time_serie)-1)] + noise)
  }

}



