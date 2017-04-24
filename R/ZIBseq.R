ZIBseq <-
function(data,outcome,transform=F,alpha=0.05) {


genodata=data     ## data is a matrix with size n*p, where n is the samoke size and p is the number of features

## must contain the column names
Y=outcome        ## Y is a vector with length n

useF=which(colSums(genodata)>2*dim(genodata)[1])    ## remove features with totle counts less than 2 times sample size

X=genodata[,useF]

ST=rowSums(X)   ## calculation the sample total

P=dim(X)[2]


##-------------------------------------------------------------------------------------


beta=matrix(data=NA,P,2)

for (i in 1:P)

{x.prop=X[,i]/ST
    
    if (transform==T)
    
    {x.prop=sqrt(x.prop) }
    
    
    bereg=gamlss(x.prop~Y,family=BEZI(sigma.link="identity"),trace=FALSE,control = gamlss.control(n.cyc = 100))
    
    out=summary(bereg)
    
    beta[i,]=out[2,c(1,4)]}  ## get all the coefficients and P values in beta regression

pvalues=beta[,2]

qvalues=calc_qvalues(pvalues)

sig=which(qvalues<alpha)

sigFeature=colnames(X)[sig]

    list(sigFeature=sigFeature,useFeature=P,qvalues=qvalues, pvalues=pvalues) }


