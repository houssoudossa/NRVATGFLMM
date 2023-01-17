\name{kin2}
\alias{kin2}
%\docType{data}
\title{
  kinship matrix
%%   ~~ data name/kind ... ~~
}
\description{
 Example of block relationship matrices of all the members in the pedigree. 
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("kin2")}
\format{
  A block matrices of 600 rows and 600 columns of NRVATGFLMM.
}
\details{
kin2 is a positive semi-definite relationship matrix (e.g. kinship matrix in genetic association studies). The rownames and colnames of this block matrices must at least include all samples as specified in the \code{ped} column of the data frame \code{Ped}.
%%  ~~ If necessary, more details than the __description__ above ~~
}
\seealso{
\code{\link{Ped}}, \code{\link{cov}}, \code{\link{geno}}, and \code{\link{snpPos}}
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
