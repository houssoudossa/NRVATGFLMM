\name{NRVAT}
\alias{NRVAT}
\title{
  Novel Rave Variants Association Test
}
\description{
Computes the Rare Variants Association Test (NRVAT) using the gaussian copula to model the relationship between relatives.
}
\usage{
NRVAT(n, gamm, heritability, a1 = 1, b1 = 25, Ped, covariates, geno, kin2)
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
}

\value{
 \code{NRVAT} returns a data frame with the following components:
 
 \item{Intercept}{ Estimated value of the Intercept.
}
  \item{gam1}{ Estimated value of the coeffeicient of the first covariate.
}
  \item{gam2}{ Estimated value of the coeffeicient of the second covariate.
}
  \item{h.S_estim}{ Estimated value of the heritability parameter.
}
  \item{NRVAT_S_obsL}{ Linear kernel test score statistic.
}
  \item{NRVAT_p.valueL}{ Linear kernel test p-value.
}
  \item{NRVAT_S_obsQ}{ Quadratic kernel test score statistic.
}
  \item{NRVAT_p.valueQ}{ Quadratic kernel test p-value.
}
  \item{NRVAT_S_obsIB}{ IBS kernel test score statistic.
}  
  \item{NRVAT_p.valueIB}{ IBS kernel test p-value.
}
  \item{NRVAT_S_obsG}{ Gaussian kernel test score statistic.
}
  \item{NRVAT_p.valueG}{ Gaussian kernel test p-value.
}
  \item{NRVAT_S_obsP}{ Polynomial kernel test score statistic.
}
  \item{NRVAT_p.valueP}{ Polynomial kernel test p-value.
}
}
\references{
}
\author{
DOSSA H. Roland G., Lakhal-Chaieb Lajmi, Oualkacha Karim 
}
\seealso{
\code{\link{CBGF_Beta_Only}}, and \code{\link{CBGF_Fixed}}
}
\examples{
data(Ped)
data(kin2)
data(geno)
data(cov)
attach(Ped)
attach(geno)
attach(cov)
members_fam=rep(c(3,4,8),40)
cov_coef = c(-2,1,1)
h.s=0.2

NRVATGFLMM_NRVAT=NRVAT(n=members_fam, gamm = cov_coef, heritability=h.s, 
  a1=1, b1=25, Ped=Ped, covariates=cov, geno=geno, kin2=kin2)

NRVATGFLMM_NRVAT
}
