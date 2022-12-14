% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_markers.R
\name{filter_snps}
\alias{filter_snps}
\title{Filter SNPs based on missingness, minor-allele frequency, and LD}
\usage{
filter_snps(
  x,
  r2.max,
  maf.min,
  indiv.miss.max,
  snp.miss.max,
  check.matrix = TRUE
)
}
\arguments{
\item{x}{An n x p marker matrix of n individuals and p markers coded as -1, 0, and 1
for homozygous alternate, heterozygous, and homozygous reference.}

\item{r2.max}{The maximum acceptable value of r2 (i.e. LD) between any two markers.}

\item{maf.min}{The minimum minor allele frequency to retain a snp.}

\item{snp.miss.max}{The maximum missingness to retain a snp}

\item{check.matrix}{Logical. Should the marker matrix 'x' be checked? Use check.matrix = FALSE
if 'x' contains imputed decimal genotypes.}

\item{entry.miss.max}{The maximum missingness to retain an individual}
}
\value{
A filtered snp matrix
}
\description{
Filter SNPs based on missingness, minor-allele frequency, and LD
}
