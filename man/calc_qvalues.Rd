\name{calc_qvalues}
\alias{calc_qvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
a function used to calculate q values
%%  ~~function to do ... ~~
}
\description{
Estimates their q-values based on a list of p-values resulting from the simultaneous testing of many hypothesis. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calc_qvalues(pvalues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalues}{
input the p value
%%     ~~Describe \code{pvalues} here~~
}
}
\details{To control the false discovery rate(FDR), q-value has been widely accepted as an alternative approach for multiple hypothesis testing correction in recent years.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
qvalues
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{http://bioconductor.org/packages/release/bioc/html/qvalue.html
%% ~put references to the literature/web site here ~
}
\author{chen hongliang
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (pvalues) 
{
    nrows = length(pvalues)
    lambdas <- seq(0, 0.95, 0.01)
    pi0_hat <- array(0, dim = c(length(lambdas)))
    for (l in 1:length(lambdas)) {
        count = 0
        for (i in 1:nrows) {
            if (pvalues[i] > lambdas[l]) {
                count = count + 1
            }
            pi0_hat[l] = count/(nrows * (1 - lambdas[l]))
        }
    }
    f <- unclass(smooth.spline(lambdas, pi0_hat, df = 3))
    f_spline <- f$y
    pi0 = f_spline[length(lambdas)]
    ordered_ps <- order(pvalues)
    pvalues <- pvalues
    qvalues <- array(0, dim = c(nrows))
    ordered_qs <- array(0, dim = c(nrows))
    ordered_qs[nrows] <- min(pvalues[ordered_ps[nrows]] * pi0, 
        1)
    for (i in (nrows - 1):1) {
        p = pvalues[ordered_ps[i]]
        new = p * nrows * pi0/i
        ordered_qs[i] <- min(new, ordered_qs[i + 1], 1)
    }
    for (i in 1:nrows) {
        qvalues[ordered_ps[i]] = ordered_qs[i]
    }
    return(qvalues)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
