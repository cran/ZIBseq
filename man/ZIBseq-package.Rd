\name{ZIBseq-package}
\alias{ZIBseq-package}
\docType{package}
\title{
Identify differentially abundant features 
}
\description{
    \packageDescription{ZIBseq}
    \packageIndices{ZIBseq}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
 \packageAuthor{ZIBseq}

 Maintainer: \packageMaintainer{ZIBseq}
}
\references{Peng Xiaoling, Li Gang, and Liu Zhenqiu. Journal of Computational Biology. January 2016, 23(2): 102-110. doi:10.1089/cmb.2015.0157.
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link{ZIBseq}} ~~
}
\examples{
\dontrun{
data(testdata)
x=testdata[,9:248]
p=dim(x)[2]
for (i in 1:p){x[,i]=as.numeric(as.character(x[,i]))}
gr=testdata[,2]
gr=as.numeric(gr)
gr[which(gr<4)]=0
gr[which(gr==4)]=1
result=ZIBseq(data=x,outcome=gr)
}
}