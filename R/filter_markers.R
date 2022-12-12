#' Calculate minor allele frequency
#'
#' @param x An n x p marker matrix of n individuals and p markers coded as -1, 0, and 1
#' for homozygous alternate, heterozygous, and homozygous reference.
#' @param check.matrix Logical. Should the marker matrix 'x' be checked? Use check.matrix = FALSE
#' if 'x' contains imputed decimal genotypes.
#'
#' @export
#'
calc_maf <- function(x, check.matrix = TRUE) {

  ## Error checking
  stopifnot(is.logical(check.matrix))
  if (check.matrix) stopifnot(check_marker_matrix(x))

  # Filter for MAF
  af <- colMeans(x + 1, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)

  return(maf)

}


#' Prune SNPs based on linkage disequilibrium
#'
#' @param x An n x p marker matrix of n individuals and p markers coded as -1, 0, and 1
#' for homozygous alternate, heterozygous, and homozygous reference.
#' @param r2max The maximum acceptable value of r2 (i.e. LD) between any two markers.
#' @param check.matrix Logical. Should the marker matrix 'x' be checked? Use check.matrix = FALSE
#' if 'x' contains imputed decimal genotypes.
#'
#' @return
#' A marker matrix without markers in high LD
#'
#' @export
#'
prune_LD <- function(x, r2.max = 0.80, check.matrix = TRUE) {

  ## Error checking
  # Check the marker matrix
  stopifnot(is.logical(check.matrix))
  if (check.matrix) stopifnot(check_marker_matrix(x))
  # check r2max
  stopifnot(r2.max >= 0 & r2.max <= 1)

  # calculate minor allele frequency
  maf <- calc_maf(x)

  # calculate the correlation between all markers; square it
  if (any(is.na(x))) {
    all_marker_r <- cor(x, use = "pairwise.complete.obs")^2

  } else {
    all_marker_r <- cor(x)^2
  }

  # Set the lower half (including diagonal) to NA
  all_marker_r1 <- all_marker_r
  all_marker_r1[lower.tri(all_marker_r1, diag = TRUE)] <- NA

  # Get a matrix of those entries that are elevated
  elevated_marker_r <- all_marker_r1 > r2.max
  # Get the coordinates of those entries
  which_elevated <- which(x = elevated_marker_r, arr.ind = TRUE)

  markers_remove <- character()
  # While loop
  i = 1
  while(nrow(which_elevated) > 0) {

    # Subset the first row
    coords <- which_elevated[1,]
    # Extract the coordinate
    r2_coord <- all_marker_r1[coords[1], coords[2], drop = FALSE]

    # marker pair
    markers <- unlist(dimnames(r2_coord))
    # Identify the marker with higher MAF
    higher_maf <- which.max(maf[markers])

    # Toss that marker
    marker_remove <- names(higher_maf)
    markers_remove[i] <- marker_remove

    # Find the row/col containing this marker
    row_remove <- col_remove <- which(row.names(all_marker_r1) == marker_remove)

    which_elevated <- subset.matrix(which_elevated, which_elevated[,"row"] != row_remove & which_elevated[,"col"] != col_remove)

    # advance i
    i <- i + 1

  }

  # Remove the markers from the marker matrix
  cols_keep <- setdiff(seq_len(ncol(x)), which(colnames(x) %in% markers_remove))
  x[,cols_keep,drop = FALSE]

}


#' Filter SNPs based on missingness, minor-allele frequency, and LD
#'
#' @param x An n x p marker matrix of n individuals and p markers coded as -1, 0, and 1
#' for homozygous alternate, heterozygous, and homozygous reference.
#' @param r2.max The maximum acceptable value of r2 (i.e. LD) between any two markers.
#' @param maf.min The minimum minor allele frequency to retain a snp.
#' @param entry.miss.max The maximum missingness to retain an individual
#' @param snp.miss.max The maximum missingness to retain a snp
#' @param check.matrix Logical. Should the marker matrix 'x' be checked? Use check.matrix = FALSE
#' if 'x' contains imputed decimal genotypes.
#'
#' @return
#' A filtered snp matrix
#'
#' @export
#'
filter_snps <- function(x, r2.max, maf.min, indiv.miss.max, snp.miss.max, check.matrix = TRUE) {

  ## Error checking
  # Check the marker matrix
  stopifnot(is.logical(check.matrix))
  if (check.matrix) stopifnot(check_marker_matrix(x))

  ## Filter on missingness
  if (!missing(snp.miss.max)) {
    stopifnot(snp.miss.max >= 0 & snp.miss.max <= 1)

    snp_missing <- colMeans(is.na(x))
    x1 <- x[,snp_missing <= snp.miss.max, drop = FALSE]

  } else {
    x1 <- x
  }

  if (!missing(indiv.miss.max)) {
    stopifnot(indiv.miss.max >= 0 & indiv.miss.max <= 1)

    indiv_missing <- rowMeans(is.na(x1))
    x2 <- x1[indiv_missing <= indiv.miss.max, , drop = FALSE]

  } else {
    x2 <- x1
  }

  ## Filter on maf
  if (!missing(maf.min)) {
    stopifnot(maf.min >= 0 & maf.min <= 1)

    # calculate minor allele frequency
    maf <- calc_maf(x2)
    x3 <- x2[,maf >= maf.min, drop = FALSE]

  } else {
    x3 <- x2
  }

  if (!missing(r2.max)) {
    stopifnot(r2.max >= 0 & r2.max <= 1)
    # Prune on LD
    x4 <- prune_LD(x = x3, r2.max = r2.max)

  } else {
    x4 <- x3
  }

  # Return x4
  return(x4)

}


