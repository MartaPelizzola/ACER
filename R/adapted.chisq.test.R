adapted.chisq.test <- function(freq, coverage, Ne, gen, poolSize=NULL, mincov=1, MeanStart=TRUE, IntGen=FALSE, TA=FALSE, RetVal=0){
  
  if(!RetVal==0 && !RetVal==1 && !RetVal==2){
    stop("(", RetVal, ") is not a valid choice for RetVal.
         RetVal needs to be 0 (p-value) or 1 (test statistic) or 2 (test statistic an p-value).")
  }
  
  if(IntGen==FALSE && length(gen)>2){
    stop("IntGen is set to FALSE: Only two generations can be considered. 
         For more than one replicate the CMH test must be used.")
  }
  
  if(IntGen==TRUE && length(unique(gen))<length(gen)){
    stop("The values of gen are not all different. 
         For more than one replicate the CMH test must be used.")
  }
  
  if(sum(is.na(freq))>0){
    stop("Allele frequency matrix has missing values.")
  }
  
  if(sum(is.na(coverage))>0){
    stop("Coverage matrix has missing values.")
  }
  
  if(sum(freq<0)>0) {
    stop("Negative allele frequencies are not allowed.")
  }
  
  if(sum(coverage<=0)>0) {
    stop("Negative values and 0 are not allowed for coverage.")
  }
  
  if(!is.null(poolSize) && sum(poolSize<=0)>0) {
    stop("Negative values and 0 are not allowed for poolSize.")
  }
  
  if(IntGen==TRUE && length(unique(gen))==2){
    warning("IntGen is set to TRUE, but only two time points are given. They will be used
            as first and last time point and no intermediate generations will be considered.")
    IntGen=FALSE
  }
  
  if (TA==TRUE && IntGen==FALSE){
    warning("Is it only possible to do Taylor approximation of the variance of the test 
            statistic if intermediate generatons are given. TA option will be ignored")
  }
  
  if(length(mincov <- as.numeric(mincov)) != 1 ) {
    stop("Length of 'mincov' (", length(mincov), ") has to be equal to '1'.")
  }
  
  if(is.na(mincov) | mincov < 1 ) {
    stop("'mincov' (", mincov, ") has to be >= 1.")
  }
  
  if (sum(coverage<mincov)>0) {
    warning("'NA' will be returned in the entries where 'coverage' is smaller than 'mincov' (", mincov, ") ")
  }
  
  if(is.vector(freq))
    freq <- matrix(freq, nrow=1)
  if(is.vector(coverage))
    coverage <- matrix(coverage, nrow=1)
  
  if(!identical(dim(freq), dim(coverage)))
    stop("The dimensions of 'freq' (", dim(freq), ") and 'coverage' (", dim(coverage), ") have to be identical.")
  
  if(!missing(Ne)){
    if (!is.integer(Ne)){
      Ne <- as.integer(Ne)
      warning("Ne value(s) which are not integer are converted to integer")
    }
  }
    
  npop <- ncol(freq)
  if(npop == 1)
    stop("Allele frequencies of at least two populations need to be provided.")
  if(npop %% length(gen) != 0) 
    stop("The number of populations (", npop, ") has to be a multiple of the length of 'gen'(", length(gen), ").")
  
  popInfo <- data.table(pop=1:ncol(freq), gen=gen)
  
  ng <- length(unique(gen))
  
  mask <- rowSums(coverage < mincov)
  stat <- NA
  
  ####### no drift
  if(length(unique(popInfo$gen)) == 1) {
    
    if(!missing(Ne))
      warning("Value of 'Ne' will be ignored because no random genetic drift is assumed.")
    
    # individual sequencing - simple chi-squared test
    if(is.null(poolSize)) {
      x1 <- coverage[,1]
      x2 <- coverage[,2]
      n <- x1 + x2
      x11 <- freq[,1] * x1
      if(sum(x11==0)>0) {
        x11[x11==0] <- rep(1,sum(x11==0))
        warning("The counts that equal 0 or equal the coverage of the considered 
                locus are changed to 1 or to coverage-1 respectively.")
      }
      if(sum(x11==x1)>0){
        x11[x11==x1] <- x1[x11==x1]-1
        warning("The counts that equal 0 or equal the coverage of the considered 
                locus are changed to 1 or to coverage-1 respectively.")
      } 
      x12 <- x1 - x11
      x21 <- freq[,2] * x2
      x22 <- x2 - x21
      
      x_1 <- x11 + x21
      x_2 <- x12 + x22
      
      stat <- (n*((x11/x1)-(x21/x2))^2)/((x_1*x_2)/(x1*x2))
      
      # pooled sequencing
      } else if(npop %% length(poolSize) == 0) {
        R1 <- coverage[,1]
        R2 <- coverage[,2]
        
        x11R <- freq[,1] * R1
        if(sum(x11R==0)>0) {
          x11R[x11R==0] <- rep(1,sum(x11R==0))
          warning("The counts that equal 0 or equal the coverage of the considered 
                  locus are changed to 1 or to coverage-1 respectively.")
        }
        if(sum(x11R==R1)>0){
          x11R[x11R==R1] <- R1[x11R==R1]-1
          warning("The counts that equal 0 or equal the coverage of the considered 
                  locus are changed to 1 or to coverage-1 respectively.")
        } 
        
        x12R <- R1 - x11R
        x21R <- freq[,2] * R2
        x22R <- R2 - x21R
        
        x1 <- poolSize[1]
        x2 <- poolSize[2]
        
        stat <- (((x11R/R1)-(x21R/R2))^2)/(((x11R*x12R)/(R1^3))*(1+(R1-1)/x1) + ((x21R*x22R)/(R2^3))*(1+(R2-1)/x2))
        } else {
          stop("The number of populations (", npop, ") has to be a multiple of the length of 'poolSize' (", length(poolSize), ").")
        }
    
    ####### with drift
      } else {
        
        # individual sequencing
        if(is.null(poolSize)) {
          ming <- min(gen)
          maxg <- max(gen)
          x1 <- coverage[,popInfo$gen == ming]
          x2 <- coverage[,popInfo$gen == maxg]
          
          x11 <- freq[,popInfo$gen == ming] * x1
          if(sum(x11==0)>0) {
            x11[x11==0] <- rep(1,sum(x11==0))
            warning("The counts that equal 0 or equal the coverage of the considered 
                    locus are changed to 1 or to coverage-1 respectively.")
          }
          if(sum(x11==x1)>0){
            x11[x11==x1] <- x1[x11==x1]-1
            warning("The counts that equal 0 or equal the coverage of the considered 
                    locus are changed to 1 or to coverage-1 respectively.")
          } 
          
          x12 <- x1 - x11
          x21 <- freq[,popInfo$gen == maxg] * x2
          x22 <- x2 - x21
          
          sigd <- x11/x1*(1-x11/x1)*(1-(1-1/(2*Ne))^(maxg-ming))
          
          Ep2 <- (x11/x1 + x21/x2)/2
          if(MeanStart==FALSE) { Ep2 <- x11/x1 }
          
          # intermediate generations
          if (IntGen == TRUE){
            
            freq_int <- freq[,!(popInfo$gen %in% c(ming, maxg))]
            coverage_int <- coverage[,!(popInfo$gen %in% c(ming, maxg))]
            x_int1 <- freq_int * coverage_int
            
            if(sum(x_int1==0)>0) {
              x_int1[x_int1==0] <- rep(1,sum(x_int1==0))
              warning('The counts of the intermediate generations assuming values 0 or equal the coverage of the considered 
                      locus are changed to 1 and to coverage-1 respectively.')
            }
            if(sum(x_int1==coverage_int)>0){
              x_int1[x_int1==coverage_int] <- coverage_int[x_int1==coverage_int]-1
              warning('The counts of the intermediate generations assuming values 0 or equal the coverage of the considered 
                      locus are changed to 1 and to coverage-1 respectively.')
            } 
            
            if (is.matrix(freq_int)==FALSE & length(freq_int)==nrow(freq)) {
                x_int1 <- matrix(x_int1, ncol=1)
                coverage_int <- matrix(coverage_int, ncol=1)
              }
            
            if (MeanStart==TRUE) {
              allCountsVec <- cbind(as.numeric(x11),x_int1,as.numeric(x21))
              allNVec <- cbind(as.numeric(x1),coverage_int,as.numeric(x2))
              allFreqVec <- allCountsVec/allNVec
              Ep2 <- rowMeans(allFreqVec)
            }
            
            if (TA==FALSE){
              if (nrow(freq)==1) {sigd <- sum((c(x11,x_int1)*(c(x1,coverage_int)-c(x11,x_int1))/c(x1,coverage_int)^2)*(1-(1-1/(2*Ne))^c(sort(gen)[2]-sort(gen)[1],sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)])))
              }
              else {sigd <- rowSums((cbind(x11,x_int1)*(cbind(x1,coverage_int)-cbind(x11,x_int1))/cbind(x1,coverage_int)^2)*(matrix(rep(1-(1-1/(2*Ne))^c(sort(gen)[2]-sort(gen)[1],sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)]),length(x1)),nrow = length(x1), byrow = TRUE)))}
            } else {
              if (nrow(freq)==1) {sigd <- sum((c(x11,x_int1)*(c(x1,coverage_int)-c(x11,x_int1))/c(x1,coverage_int)^2)*(c(sort(gen)[2]-sort(gen)[1],sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)]))/((2*Ne)^1))}
              else {sigd <- rowSums((cbind(x11,x_int1)*(cbind(x1,coverage_int)-cbind(x11,x_int1))/cbind(x1,coverage_int)^2)*(matrix(rep(c(sort(gen)[2]-sort(gen)[1],sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)]),length(x1)),nrow = length(x1), byrow = TRUE))/((2*Ne)^1)) }
            }
            }
          
          stat <- ((x11*x22-x12*x21)^2)/(x2^2*x11*x12/x1+x1^2*x2*(Ep2*(1-Ep2)+(x2-1)*sigd))
          
          # pooled sequencing
          } else if(npop %% length(poolSize) == 0) {
            ming <- min(gen)
            maxg <- max(gen)
            R1 <- coverage[,popInfo$gen == ming]
            R2 <- coverage[,popInfo$gen == maxg]
            
            x11R <- freq[,popInfo$gen == ming] * R1
            if(sum(x11R==0)>0) {
              x11R[x11R==0] <- rep(1,sum(x11R==0))
              warning("The counts that equal 0 or equal the coverage of the considered 
                      locus are changed to 1 or to coverage-1 respectively.")
            }
            if(sum(x11R==R1)>0){
              x11R[x11R==R1] <- R1[x11R==R1]-1
              warning("The counts that equal 0 or equal the coverage of the considered 
                      locus are changed to 1 or to coverage-1 respectively.")
            } 
            
            x12R <- R1 - x11R
            x21R <- freq[,popInfo$gen == maxg] * R2
            x22R <- R2 - x21R
            
            x1 <- poolSize[1]
            x2 <- poolSize[2]
            
            Ep2 <- (x11R/R1+x21R/R2)/2
            if(MeanStart==FALSE){Ep2 <- x11R/R1}
            
            sigd <- (x11R/R1)*(1-(x11R/R1))*(1-(1-1/(2*Ne))^(maxg-ming))
            
            # intermediate generations
            if (IntGen == TRUE){
              
              freq_int <- freq[,!(popInfo$gen %in% c(ming, maxg))]
              coverage_int <- coverage[,!(popInfo$gen %in% c(ming, maxg))]
              
              x_int1 <- freq_int * coverage_int
              
              if(sum(x_int1==0)>0) {
                x_int1[x_int1==0] <- rep(1,sum(x_int1==0))
                warning('The counts of the intermediate generations assuming values 0 or equal the coverage of the considered 
                        locus are changed to 1 and to coverage-1 respectively.')
              }
              if(sum(x_int1==coverage_int)>0){
                x_int1[x_int1==coverage_int] <- coverage_int[x_int1==coverage_int]-1
                warning('The counts of the intermediate generations assuming values 0 or equal the coverage of the considered 
                        locus are changed to 1 and to coverage-1 respectively.')
              } 
              
              if (is.matrix(freq_int)==FALSE & length(freq_int)==nrow(freq)) {
                x_int1 <- matrix(x_int1, ncol=1)
                coverage_int <- matrix(coverage_int, ncol=1)
              }
              
              if (MeanStart==TRUE) {
                allCountsVec <- cbind(as.numeric(x11R),x_int1,as.numeric(x21R))
                allNVec <- cbind(as.numeric(R1),coverage_int,as.numeric(R2))
                allFreqVec <- allCountsVec/allNVec
                Ep2 <- rowMeans(allFreqVec)
              }
              
              if (TA==FALSE){
                if (nrow(freq)==1) {sigd <- sum((c(x11R,x_int1)*(c(R1,coverage_int)-c(x11R,x_int1))/c(R1,coverage_int)^2)*(1-(1-1/(2*Ne))^c(sort(gen)[2]-sort(gen)[1],sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)]))) }
                else { sigd <- rowSums((cbind(x11R,x_int1)*(cbind(R1,coverage_int)-cbind(x11R,x_int1))/cbind(R1,coverage_int)^2)*(matrix(rep(1-(1-1/(2*Ne))^c(sort(gen)[2]-sort(gen)[1],sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)]),length(R1)),nrow = length(R1), byrow = TRUE))) }
              } else {
                if (nrow(freq)==1) {sigd <- sum((c(x11R,x_int1)*(c(R1,coverage_int)-c(x11R,x_int1))/c(R1,coverage_int)^2)*(c(sort(gen)[2]-sort(gen)[1],sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)]))/((2*Ne)^1)) }
                else { sigd <- rowSums((cbind(x11R,x_int1)*(cbind(R1,coverage_int)-cbind(x11R,x_int1))/cbind(R1,coverage_int)^2)*(matrix(rep(c(sort(gen)[2]-sort(gen)[1],sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)]),length(R1)),nrow = length(R1), byrow = TRUE))/((2*Ne)^1)) }
              }
              }
            
            stat <- (x11R*x22R-x12R*x21R)^2/(R2^2*(R1*(x11R/R1)*(1-(x11R/R1))*(1+(R1-1)/x1))+R1^2*(R2*Ep2-R2*Ep2^2 + R2*(R2-1)/x2*(Ep2*(1-Ep2)+(x2-1)*sigd)))
            
            } else {
              stop("The number of populations (", npop, ") has to be a multiple of the length of 'poolSize' (", length(poolSize), ").")
            }
      }
  
  # convert summary statistic to p-value and return it
  stat[mask>0] <- NA
  res <- pchisq(stat, df=1, lower.tail=FALSE) 
  if(!RetVal==0){
    if(RetVal==1){
      res <- stat 
    }else if(RetVal==2){
      res <- cbind(stat,pchisq(stat, df=1, lower.tail=FALSE))
      colnames(res) <- c("test_statistic", "p.value")
    }
  }
  
  return(res)
  }
