#' Hidden functions
#'
check_marker_matrix <- function(x) {
  # First check if it's a matrix
  if (!is.matrix(x)) {
    return(FALSE)

  } else if (!all(na.omit(unique(as.vector(x))) %in% c(-1, 0, 1), na.rm = TRUE)) {
    return(FALSE)

  } else if (any(sapply(dimnames(x), is.null))) {
    return(FALSE)

  } else {
    return(TRUE)

  }

}
