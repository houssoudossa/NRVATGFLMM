\name{geno}
\alias{geno}

\title{
  geno
}
\description{
Example genotype data frame.
}
\usage{data("geno")}
\format{
  A data frame with 600 rows and 102 columns of NRVATGFLMM.
}
\details{
A data frame containing the genotype information. It is a matrix with genotypes for subjects
(rows) at each variant position (columns). The first two columns are needed to be named “ped” and “person”, which are used to match subjects to their data in the pedigree data.frame. The genotypes
are coded as 0, 1, 2 for autosomal markers (generally a count of the number of the minor alleles).
}
\seealso{
\code{\link{Ped}}, \code{\link{cov}}, \code{\link{snpPos}}, and \code{\link{Kin2}}
}
