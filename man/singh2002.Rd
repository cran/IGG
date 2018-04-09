\name{singh2002}
\alias{singh2002}

\title{Prostate Cancer Study of Singh et al. (2002)}
\description{ Gene expression dat (6033 genes for 102 samples) from the microarray study of Singh et al. (2002). Also available in the \code{sda} package.
 }
\usage{data(singh2002)}
\format{
A list of two components.
\describe{
    \item{\code{x}:}{is a 102 \eqn{\times} 6033 matrix containing the expression levels. The rows contain the samples and the columns the genes.}
    \item{\code{y}:}{is a factor containing the diagnosis for each sample ("cancer" or "healthy").}
  }
}

\details{This data set contains measurements of the gene expression of 6033 genes for 102 observations. The first 52 rows are for the cancer patients and the last 50 rows are for the normal control subjects.
}
\source{ 
The data is described in Singh et al. (2001) and are provided in exactly the form as used by Efron (2010).

}
\references{
Singh, D., Febbo, P.G., Ross, K., Jackson, D.G., Maonla, J., Ladd, C., Tamayo, P., Renshaw, A.A., D'Amico, A.V., Richie, J.P., Lander, E.S., Loda, M., Kantoff, P.W., Golub, T.R., and Sellers, W.R. (2002). "Gene expression correlates of clinical prostate cancer behavior." \emph{Cancer Cell}, \bold{1}(2): 203-209.

Efron, B. (2010). "The Future of Indirect Evidence." \emph{Statistical Science}, \bold{25}(2): 145-157.
}