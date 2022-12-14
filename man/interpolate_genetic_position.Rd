% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpolate_genetic_positions.R
\name{interpolate_genetic_position}
\alias{interpolate_genetic_position}
\title{Predict marker genetic position using linear interpolation}
\usage{
interpolate_genetic_position(map, cMMb)
}
\arguments{
\item{map}{A data frame with three columns: marker name, chromosome, and physical
position. These are markers whose genetic position is unknown}

\item{cMMb}{A data frame with six columns: chromosome, left flank physical position,
right flank physical position, local recombination rate (in cM/Mb) in that flank,
left flank genetic position, and right flank genetic position.}
}
\description{
Uses linear interpolation and local recombination rate information to predict
the unknown genetic position of markers based on their known physical positions.
}
