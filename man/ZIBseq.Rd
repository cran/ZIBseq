\name{ZIBseq}
\alias{ZIBseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conducts the zero-inflated beta regression based on the general count \code{data} and categorical vector \code{outcome}.
%%  ~~function to do ... ~~
}
\description{
zero-inflated beta regression
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ZIBseq(data, outcome, transform = F, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a matrix records the count data
%%     ~~Describe \code{data} here~~
}
  \item{outcome}{
a categorical vector of a specific kind of clinical condition
%%     ~~Describe \code{outcome} here~~
}
  \item{transform}{
square-root transform of the compositional matrix
%%     ~~Describe \code{transform} here~~
}
  \item{alpha}{
customized threshold while calculating q values
%%     ~~Describe \code{alpha} here~~
}
}
\details{The function takes the sparse nature of metagenomics data into account and handle the compositional data efficiently.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
    \item{sigFeature}{output the significant feature}
    \item{useFeature}{features being concerned}
    \item{qvalue}{qvalue}
    \item{pvalue}{pvalue}
%% ...
}
\references{Peng Xiaoling, Li Gang, and Liu Zhenqiu. Journal of Computational Biology. January 2016, 23(2): 102-110. doi:10.1089/cmb.2015.0157.
%% ~put references to the literature/web site here ~
}
\author{Hongliang Chen
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{calc_qvalues}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
