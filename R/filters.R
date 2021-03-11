

#' Moving average filter.
#'
#' This function return the data smoothed using the basic moving average
#' algorithm. For each chunk of data of size equal to the buffer_size parameter
#' is calculated the average and this value is used as the i term of the newly
#' smoothed data.
#' zero padding is applied for initial and final values
#'
#' @param raw_data Data upon which the algorithm is applied
#' @param buffer_size number of points the algorithm use to compute the average
#'
#' @return Smoothed data using moving average algorithm
#'
#' @examples
#' raw_data = c(1:100)
#' smoothed_data = MA(raw_data)
#'
#' @importFrom stats median
#'
#' @export
#'

MA <- function(raw_data, buffer_size = 5) {
  if (buffer_size %% 2 == 0) {
    return(warning("odd buffer size is required in order to compute simmetrical corrections"))
  } else if (buffer_size == 1){
    return(warning("three point filter is the minimum value"))
  } else if (!(is.numeric(as.integer(buffer_size))) | is.na(as.integer(buffer_size))){
    return(warning("buffer size should be an odd integer"))
  } else if (typeof(raw_data) == "list"){
    raw_data = raw_data[[1]]
  }
  smoothed_data <- raw_data

  for (i in 1:length(raw_data)) {
    if (i < median(1:buffer_size)) {
      smoothed_data[i] <- sum(raw_data[1:(i + floor(buffer_size/2))])/buffer_size
    } else if (i > (length(raw_data) - median(1:buffer_size) + 1)) {
      smoothed_data[i] <- sum(raw_data[(i - floor(buffer_size/2)):length(raw_data)])/buffer_size
    } else {
      smoothed_data[i] <- sum(raw_data[(i - floor(buffer_size/2)):(i + floor(buffer_size/2))])/buffer_size
    }
  }
  return(smoothed_data)
}



#' Running median smoothing.
#'
#' This function return the data smoothed using the running median
#' algorithm. For each chunk of data of size equal to the buffer_size parameter
#' is calculated the median and this value is used as the i term of the newly
#' smoothed data.
#' For initial and final values zero padding is applied.
#'
#'
#' @param raw_data Data upon which the algorithm is applied
#' @param buffer_size number of points the algorithm use to compute the median
#'
#' @return Smoothed data using running median algorithm
#'
#' @examples
#' raw_data = c(1:100)
#' smoothed_data = RMS(raw_data)
#'
#' @importFrom stats median
#'
#' @export
#'

RMS <- function(raw_data, buffer_size = 5) {
  if (buffer_size %% 2 == 0) {
    return(warning("odd buffer size is required in order to compute simmetrical corrections"))
  } else if (buffer_size == 1){
    return(warning("three point filter is the minimum value"))
  } else if (!(is.numeric(as.integer(buffer_size))) | is.na(as.integer(buffer_size))){
    return(warning("buffer size should be an odd integer"))
  } else if (typeof(raw_data) == "list"){
    raw_data = raw_data[[1]]
  }
  smoothed_data <- raw_data

  for (i in 1:length(raw_data)) {
    if (i < median(1:buffer_size)) {
      median_buffer <- raw_data[1:(i + floor(buffer_size/2))]
      if ((length(median_buffer) %% 2 == 0)) {
        median_buffer <- append(median_buffer, 0)
      }
      smoothed_data[i] = median(median_buffer)
    } else if (i > (length(raw_data) - median(1:buffer_size) + 1)) {
      median_buffer <- raw_data[(i - floor(buffer_size/2)):length(raw_data)]
      if ((length(median_buffer) %% 2 == 0)) {
        median_buffer <- append(median_buffer, 0)
      }
      smoothed_data[i] <- median(median_buffer)
    } else {
      smoothed_data[i] <- median(raw_data[(i - floor(buffer_size/2)):(i + floor(buffer_size/2))])
    }
  }
  return(smoothed_data)
}


#' Hann window filter.
#'
#' This function return the data smoothed using the an Hann window filter.
#' Data are smoothed using a cosine window.
#'
#' @param raw_data Data upon which the algorithm is applied
#' @param buffer_size number of points the algorithm use to compute the coefficients of the Hann window
#'
#' @return Smoothed data using Hann Window filter
#'
#' @examples
#' raw_data = c(1:100)
#' smoothed_data = Hann(raw_data)
#'
#' @importFrom stats median
#'
#' @export
#'

Hann <- function(raw_data, buffer_size = 5) {
  if (buffer_size %% 2 == 0) {
    return(warning("odd buffer size is required in order to compute simmetrical corrections"))
  } else if (buffer_size == 1){
    return(warning("three point filter is the minimum value"))
  } else if (!(is.numeric(as.integer(buffer_size))) | is.na(as.integer(buffer_size))){
    return(warning("buffer size should be an odd integer"))
  } else if (typeof(raw_data) == "list"){
    raw_data = raw_data[[1]]
  }
  smoothed_data <- raw_data

  hanncoeff <- 0.5*(1 - cos(2*pi*c(0:(buffer_size-1))/(buffer_size-1)))
  hanncoeff <- hanncoeff/sum(hanncoeff)


  for (i in 1:length(raw_data)) {
    if (i < median(1:buffer_size)) {
      data_buffer <- rev(raw_data[1:(i + floor(buffer_size/2))])
      while ((length(data_buffer) != buffer_size)) {
        data_buffer <- append(data_buffer, 0)
      }
      smoothed_data[i] <- sum(hanncoeff*rev(data_buffer))
    } else if (i > (length(raw_data) - median(1:buffer_size) + 1)) {
      data_buffer <- raw_data[(i - floor(buffer_size/2)):length(raw_data)]
      while ((length(data_buffer) != buffer_size)) {
        data_buffer <- append(data_buffer, 0)
      }
      smoothed_data[i] <- sum(hanncoeff*data_buffer)
    } else {
      smoothed_data[i] <- sum(hanncoeff*raw_data[(i - floor(buffer_size/2)):(i + floor(buffer_size/2))])
    }
  }
  return(smoothed_data)
}

#' Hamming window filter.
#'
#' This function return the data smoothed using the an Hamming window filter.
#' Data are smoothed using a cosine window with particular coefficients.
#'
#' @param raw_data Data upon which the algorithm is applied
#' @param buffer_size number of points the algorithm use to compute the coefficients of the Hann window
#'
#' @return Smoothed data using Hann Window filter
#'
#' @examples
#' raw_data = c(1:100)
#' smoothed_data = Hamming(raw_data)
#'
#' @importFrom stats median
#'
#' @export
#'

Hamming <- function(raw_data, buffer_size = 5) {
  if (buffer_size %% 2 == 0) {
    return(warning("odd buffer size is required in order to compute simmetrical corrections"))
  } else if (buffer_size == 1){
    return(warning("three point filter is the minimum value"))
  } else if (!(is.numeric(as.integer(buffer_size))) | is.na(as.integer(buffer_size))){
    return(warning("buffer size should be an odd integer"))
  } else if (typeof(raw_data) == "list"){
    raw_data = raw_data[[1]]
  }
  smoothed_data <- raw_data

  hammcoeff <- (0.54 - 0.46*cos(2*pi*c(0:(buffer_size-1))/(buffer_size-1)))
  hammcoeff <- hammcoeff/sum(hammcoeff)


  for (i in 1:length(raw_data)) {
    if (i < median(1:buffer_size)) {
      data_buffer <- rev(raw_data[1:(i + floor(buffer_size/2))])
      while ((length(data_buffer) != buffer_size)) {
        data_buffer <- append(data_buffer, 0)
      }
      smoothed_data[i] <- sum(hammcoeff*rev(data_buffer))
    } else if (i > (length(raw_data) - median(1:buffer_size) + 1)) {
      data_buffer <- raw_data[(i - floor(buffer_size/2)):length(raw_data)]
      while ((length(data_buffer) != buffer_size)) {
        data_buffer <- append(data_buffer, 0)
      }
      smoothed_data[i] <- sum(hammcoeff*data_buffer)
    } else {
      smoothed_data[i] <- sum(hammcoeff*raw_data[(i - floor(buffer_size/2)):(i + floor(buffer_size/2))])
    }
  }
  return(smoothed_data)
}


#' Sine window filter.
#'
#' This function return the data smoothed using the a sine window filter.
#'
#' @param raw_data Data upon which the algorithm is applied
#' @param buffer_size number of points the algorithm use to compute the coefficients of the Hann window
#'
#' @return Smoothed data using Hann Window filter
#'
#' @examples
#' raw_data = c(1:100)
#' smoothed_data = sine(raw_data)
#'
#' @importFrom stats median
#'
#' @export
#'

sine <- function(raw_data, buffer_size = 5) {
  if (buffer_size %% 2 == 0) {
    return(warning("odd buffer size is required in order to compute simmetrical corrections"))
  } else if (buffer_size == 1){
    return(warning("three point filter is the minimum value"))
  } else if (!(is.numeric(as.integer(buffer_size))) | is.na(as.integer(buffer_size))){
    return(warning("buffer size should be an odd integer"))
  } else if (typeof(raw_data) == "list"){
    raw_data = raw_data[[1]]
  }

  smoothed_data <- raw_data

  sinecoeff <- sin(pi*c(0:(buffer_size-1))/(buffer_size-1))
  sinecoeff <- sinecoeff/sum(sinecoeff)


  for (i in 1:length(raw_data)) {
    if (i < median(1:buffer_size)) {
      data_buffer <- rev(raw_data[1:(i + floor(buffer_size/2))])
      while ((length(data_buffer) != buffer_size)) {
        data_buffer <- append(data_buffer, 0)
      }
      smoothed_data[i] <- sum(sinecoeff*rev(data_buffer))
    } else if (i > (length(raw_data) - median(1:buffer_size) + 1)) {
      data_buffer <- raw_data[(i - floor(buffer_size/2)):length(raw_data)]
      while ((length(data_buffer) != buffer_size)) {
        data_buffer <- append(data_buffer, 0)
      }
      smoothed_data[i] <- sum(sinecoeff*data_buffer)
    } else {
      smoothed_data[i] <- sum(sinecoeff*raw_data[(i - floor(buffer_size/2)):(i + floor(buffer_size/2))])
    }
  }
  return(smoothed_data)
}
