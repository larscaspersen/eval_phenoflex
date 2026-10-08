# Normalize only characteristic temperature parameters, including indexed names.
# Keep the native kernels' legacy Celsius offset; Kelvin inputs are unchanged.
.characteristic_temperature_names <- function(names) {
  grepl("^(theta_star|theta_c)[0-9]*$", names)
}

.normalize_characteristic_temperatures <- function(parameters) {
  # Malformed inputs are left for the calling interface's existing validation.
  if (!is.numeric(parameters) || !is.null(dim(parameters)) ||
      is.null(names(parameters)) || anyNA(names(parameters)) || anyDuplicated(names(parameters)))
    return(parameters)
  celsius <- .characteristic_temperature_names(names(parameters)) &
    is.finite(parameters) & parameters >= 0 & parameters <= 20
  if (any(celsius)) parameters[celsius] <- parameters[celsius] + 273
  parameters
}
