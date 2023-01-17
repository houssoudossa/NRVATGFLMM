\name{Ped}
\alias{Ped}
\title{
  Ped
}
\description{
 Example pedigree data frame.
}
\usage{data("Ped")}
\format{
  A data frame with 600 rows and 7 columns of NRVATGFLMM.
}
\details{
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
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\seealso{
\code{\link{cov}}, \code{\link{geno}}, \code{\link{snpPos}}, and \code{\link{Kin2}}
}
