\name{SNPs.CommonSingle}
\alias{SNPs.141CommonSingle}
\alias{SNPs.142CommonSingle}
\alias{SNPs.144CommonSingle}
\alias{SNPs.146CommonSingle}
\alias{SNPs.147CommonSingle}
\alias{SNPs.150CommonSingle}
\alias{SNPs.151CommonSingle}
\docType{data}
\title{
SNP annotation from various versions of dbSNP as represented on UCSC
Common SNP table.
}
\description{
SNP annotation from various versions of dbSNP as represented on UCSC
Common SNP table. Overlap is based on genomic mappings from the
annotation package.
}
\format{
  An object of class \code{DataFrame}.  Rownames are CpG identifiers.
  There are 6 columns \code{Probe_rs}, \code{Probe_maf}, \code{CpG_rs},
  \code{CpG_maf}, \code{SBE_rs}, and \code{SBE_maf}.  \sQuote{Probe}
  indicates a SNP in the probe, \sQuote{CpG} a SNP in the CpG site and
  \sQuote{SBE} in the single base extension site.  The \code{_rs} gives
  the SNP RS identifier and the \code{_maf} gives the minor allele frequency.
  }
\details{
  In addition to the SNP information provided by Illumina, we have added
  independent information on the overlap of the EPIC (v2.0) array with various
  versions of dbSNP.  The overlap is based on the
  mappings of the array to the hg38 genome provided by Illumina.  As
  dbSNP we have used the \sQuote{Common} table from UCSC
  (ie. \sQuote{snp151Common}).  This track contains variants from dbSNP
  which have a minor allele frequency (MAF) of greater than 1 percent
  (specifically, this requires dbSNP to actually contain MAF
  information).  Furthermore, we only kept variants marked as
  \sQuote{single} (ie. standard single nucleotide changes, but not
  insertions or deletions).  Scripts for retrieving the UCSC dbSNP table
  and doing the overlap are contained in the \code{scripts} directory.
  The variants are described in 6 different columns. \code{Probe_rs}
  tells us the RS number (SNP ID number) for a SNP overlapping the
  probe, and \code{Probe_maf} is the minor allelle frequency for the SNP
  (in case multiple SNPs overlap, only one is recorded).  Similarly,
  \code{CpG_rs} describe SNPs overlapping the CpG site and \code{SBE_rs}
  the single base extension of the measured methylation loci. 
}
\source{
  UCSC Common SNP Table.
}
\examples{
SNPs.141CommonSingle
SNPs.142CommonSingle
SNPs.144CommonSingle
SNPs.146CommonSingle
SNPs.147CommonSingle
SNPs.150CommonSingle
SNPs.151CommonSingle
}
\keyword{datasets}
