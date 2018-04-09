\name{diabetes}
\alias{diabetes}

\title{Blood and other measurements in diabetics}
\description{ The \code{diabetes} data frame has 442 rows and 3 columns. These are the data used in the Efron et al. "Least Angle Regression" paper and is also available in the \code{lars} package.
 }
\usage{data(diabetes)}
\format{
This data frame consists of the following columns:
\describe{
    \item{\code{x}:}{is a design matrix with 10 columns (no interactions).}
    \item{\code{y}:}{is a numeric vector.}
    \item{\code{x2}:}{is a design matrix with 64 columns (includes interactions.}
  }
}

\details{The x matrix has been standardized to have unit L2 norm in each column and zero mean. The matrix x2 consists of x plus certain interactions.
}
\source{ 
https://cran.r-project.org/web/packages/lars/
}
\references{
Efron, B., Hastie, T., Johnstone, I., and Tibshirani, R. (2003). "Least Angle Regression" (with discussion). \emph{Annals of Statistics}, \bold{32}(2): 407-499.
}