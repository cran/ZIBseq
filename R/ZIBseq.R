require(gamlss)
library(phyloseq)
library(metagenomeSeq)

ZIBseq <-
  function(data,meta,Y,Z,transform=T,alpha=0.05) {
    
    # data: Dataframe, Sample by feature table
    # meta: Dataframe, Meta data or clinical variables for the samples including label information
    # Y: String, label name
    # Z: Srting (or string vector for multiple variables), clinical variable(s) to be used
    # transformation: Boolean, true for square root transformation
    # alpha: Float, significance level 
    
    ## Data processing

    # Handling numerical labels
    meta[,Y] = as.factor(meta[,Y])
    
    # Relative abundance transformation for beta distribution
    RA.data = t(apply(data,1,function(X){
      X/sum(X)
    }))
    
    # Square root transformation
    if(transform==T){
      RA.data = sqrt(RA.data)
    }
    
    # Dataframe for gamlss function
    df = data.frame(cbind(RA.data), meta)
    
    # Feature names for gamlss function
    Xs = colnames(RA.data)
    
    ## Zero-inflated beta regression
    p_vec = c() # Vector with p-values for each feature
    for(i in Xs){
      equat = as.formula(paste0(i,"~",paste0(c(Y,Z), collapse = " + ")))
      p.value <- tryCatch({
        bereg=gamlss(formula = equat,family=BEZI(sigma.link="identity"),trace=FALSE,control = gamlss.control(n.cyc = 100), data = df)
        out=summary(bereg)
        out[2,4] # p-value
      }, warning=function(w) {
        message("handling warning: ", conditionMessage(w))
        NA
      })
                
      p_vec <- c(p_vec, p.value)
    }
    return(p_vec)
  }

ZIBseq.multinomial <-
  function(data,meta,Y,Z,transform=T,alpha=0.05) {
    
    # data: Dataframe, Sample by feature table
    # meta: Dataframe, Meta data or clinical variables for the samples including label information
    # Y: String, label name
    # Z: Srting (or string vector for multiple variables), clinical variable(s) to be used
    # transformation: Boolean, true for square root transformation
    # alpha: Float, significance level 
    
    ## Data processing  
    # Handling numerical labels
    meta[,Y] = as.factor(meta[,Y])
    
    # Relative abundance transformation for beta distribution
    RA.data = t(apply(data,1,function(X){
      X/sum(X)
    }))
    
    # Square root transformation
    if(transform==T){
      RA.data = sqrt(RA.data)
    }
    
    # Dataframe for gamlss function
    df = data.frame(cbind(RA.data), meta)
    
    # Feature names for gamlss function
    Xs = colnames(RA.data)
    
   
    ## Zero-inflated beta regression + Likelihood ratio test
    p_vec = c() # Vector with p-values for each feature
    for(i in Xs){
      full.equation = paste0(i, " ~ ", paste0(c(Y,Z), collapse = " + "))
      reduced.equation = paste0(i, " ~ ", paste0(c(Z), collapse = " + "))
      
      full.model = gamlss(formula = as.formula(full.equation), family=BEZI(sigma.link="identity"), trace=FALSE, control = gamlss.control(n.cyc = 100), data = df)
      
      reduced.model = gamlss(formula = as.formula(reduced.equation), family=BEZI(sigma.link="identity"), trace=FALSE, control = gamlss.control(n.cyc = 100), data = df)
      
      if(deviance(reduced.model) <= deviance(full.model)){
        p.value = NA
      }else{
        p.value = LR.test(reduced.model, full.model, print = F)$p.val  
      }
      
      p_vec <- c(p_vec, p.value)
    }
    return(p_vec)
  }


