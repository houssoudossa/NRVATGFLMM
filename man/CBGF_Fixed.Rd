\name{CBGF_Fixed}
\alias{CBGF_Fixed}
\title{
  CBGF_Fixed
}
\description{
  Computes the NRVATGFLMM statistics when smoothing the genetic effect function and the genetic variant function.
}
\usage{
CBGF_Fixed(n, gamm, heritability, a1 = 1, b1 = 25, Ped, 
  covariates, geno, kin2, pos, order, beta_basis_F_BS, 
  geno_basis_F_BS, base_F_BS = "bspline", beta_basis_F_FS, 
  geno_basis_F_FS, base_F_FS = "fspline")
}
\arguments{
  \item{n}{
  A vector containing the number of members in each pedigree.
}
  \item{gamm}{
  A vector of coefficients of the covariates inclunding the intercept.
}
  \item{heritability}{
  A value (betwen 0 and 1) of the heritability parameter. It defined as the relative importance of influence genetics in phenotypic variance.
}
  \item{a1}{
  The first value of the numeric vector of length 2 defining the beta probability density function parameters on the minor allele frequencies.
}
  \item{b1}{
  The second value of the numeric vector of length 2 defining the beta probability density function parameters on the minor allele frequencies.
}
  \item{Ped}{
  A data frame containing the pedigree information with the following columns:
  
  \code{ID} { person ID
    }
    \code{ped} { pedigree ID, character or numeric allowed.
    }
    \code{person} { person ID, a unique ID within each pedigree, numeric or character allowed.
    }
    \code{father} { father ID, 0 if no father.
    }
    \code{mother} { mother ID, 0 if no mother.
    }
    \code{sex} { sex, coded as 1 for male, 2 for female.
    }
    \code{trait} { trait phenotype, either case-control status coded as 1 for affected and 0 for unaffected. Subjects with missing (NA) will be removed from the analysis.
    }
}
  \item{covariates}{
  A data frame containing the covariates information. The first two columns are needed to be named “ped” and “person”, which are used to match subjects to their data in the pedigree data frame.
}
  \item{geno}{
  A data frame containing the genotype information. It is a matrix with genotypes for subjects
(rows) at each variant position (columns). The first two columns are needed to be named “ped” and “person”, which are used to match subjects to their data in the pedigree data.frame. The genotypes are coded as 0, 1, 2 for autosomal markers (generally a count of the number of the minor alleles).

}
  \item{kin2}{
a known positive semi-definite relationship matrix (e.g. kinship matrix in genetic association studies). The rownames and colnames of these matrices must at least include all samples as specified in the ped column of the data frame \code{\link{Ped}}.
}
  \item{pos}{
  Position of the SNPs in base pairs.
}
  \item{order}{
  The order used to generate the B-spline basis.
}
  \item{beta_basis_F_BS}{ The number of B-spline basis functions used to estimate the genetic effect function.
}
  \item{geno_basis_F_BS}{ The number of B-spline basis functions used to estimate the genetic variant functions.
}
  \item{base_F_BS}{
  It is "bspline".
}
  \item{beta_basis_F_FS}{ The number of fourier basis functions used to estimate the genetic effect function
}
  \item{geno_basis_F_FS}{ The number of fourier basis functions used to estimate the genetic variant functions.
}
  \item{base_F_FS}{
  It is "fspline".
}
}

\value{
\code{CBGF_Fixed} returns a data frame with the following components:

  \item{Intercept}{ Estimated value of the Intercept.
}
  \item{gam1}{ Estimated value of the coeffeicient of the first covariate.
}
  \item{gam2}{ Estimated value of the coeffeicient of the second covariate.
}
  \item{h.S_estim}{ Estimated value of the heritability parameter.
}
  \item{StatL.Fixed_Bs}{ Linear kernel test score statistic under the B-spline basis functions.
}
  \item{p.v_obsbL.Fixed_Bs}{ Linear kernel test p-value under the B-spline basis functions.
}
  \item{StatQ.Fixed_Bs}{ Quadratic kernel test score statistic under the B-spline basis 
functions.
}
  \item{p.v_obsbQ.Fixed_Bs}{ Quadratic kernel test p-value under the B-spline basis functions.
}
  \item{StatG.Fixed_Bs}{ Gaussian kernel test score statistic under the B-spline basis functions.
}
  \item{p.v_obsbG.Fixed_Bs}{ Gaussian kernel test p-value under the B-spline basis functions.
}
  \item{StatP.Fixed_Bs}{ Polynomial kernel test score statistic under the B-spline basis functions.
}
  \item{p.v_obsbP.Fixed_Bs}{ Polynomial kernel test p-value under the B-spline basis functions.
}
  \item{StatL.Fixed_Fs}{ Linear kernel test score statistic under the fourier basis functions.
}
  \item{p.v_obsbL.Fixed_Fs}{ Linear kernel test p-value under the fourier basis functions.
}
  \item{StatQ.Fixed_Fs}{ Quadratic kernel test score statistic under the fourier basis functions.
}
  \item{p.v_obsbQ.Fixed_Fs}{ Quadratic kernel test p-value under the fourier basis functions.
}
 \item{StatG.Fixed_Fs}{ Gaussian kernel test score statistic under the fourier basis functions.
}
  \item{p.v_obsbG.Fixed_Fs}{ Gaussian kernel test p-value under the fourier basis functions.
}
  \item{StatP.Fixed_Fs}{ Polynomial kernel test score statistic under the fourier basis 
functions.
}
  \item{p.v_obsbP.Fixed_Fs}{ Polynomial kernel test p-value under the fourier basis functions.
}
}
\references{

}
\author{
  DOSSA H. Roland G., Lakhal-Chaieb Lajmi, Oualkacha Karim
}
\seealso{
  \code{\link{NRVAT}}, and \code{\link{CBGF_Beta_Only}}
}
\examples{
data(Ped)
data(kin2)
data(geno)
data(cov)
data(snpPos)
attach(Ped)
attach(geno)
attach(cov)
attach(snpPos)
members_fam=rep(c(3,4,8),40)
cov_coef = c(-2,1,1)
h.s=0.2
beta_basis_F_BS = geno_basis_F_BS = 10
beta_basis_F_FS = geno_basis_F_FS = 11
order = 4

NRVATGFLMM_Fixed=CBGF_Fixed(n=members_fam, gamm = cov_coef, 
  heritability=h.s, a1=1, b1=25, Ped=Ped, covariates=cov, 
  geno=geno, kin2=kin2, pos = snpPos, order=order, 
  beta_basis_F_BS = beta_basis_F_BS, geno_basis_F_BS = geno_basis_F_BS, 
  base_F_BS = "bspline", beta_basis_F_FS = beta_basis_F_FS, 
  geno_basis_F_FS = geno_basis_F_FS, base_F_FS ="fspline")

NRVATGFLMM_Fixed
}