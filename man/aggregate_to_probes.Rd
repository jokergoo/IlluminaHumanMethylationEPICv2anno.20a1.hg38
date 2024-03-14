\name{aggregate_to_probes}
\alias{aggregate_to_probes}
\title{
Aggregate to the probe-level
}
\description{
Aggregate to the probe-level
}
\usage{
aggregate_to_probes(x)
}
\arguments{

  \item{x}{A data frame or a matrix.}

}
\details{
If \code{x} is a data frame, it should be one of the following objects: Islands.UCSC, Locations, Other, SNPs.Illumina, SNPs.141CommonSingle, SNPs.142CommonSingle, SNPs.144CommonSingle, 
 SNPs.146CommonSingle, SNPs.147CommonSingle, SNPs.150CommonSingle, SNPs.151CommonSingle.

If \code{x} is a matrix, it should come from the analysis with the **minfi** package. If multiple rows correspond to the same
probe, the average value of rows is simply used.
}
\examples{
# There is no example
NULL
}
