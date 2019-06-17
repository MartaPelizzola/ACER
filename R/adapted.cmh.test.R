adapted.cmh.test <- function(freq, coverage, Ne, gen, repl, poolSize=NULL, mincov=1, MeanStart=TRUE, IntGen=FALSE, TA=FALSE, order=0, RetVal=0){

  if(!RetVal==0 && !RetVal==1 && !RetVal==2){
    stop("(", RetVal, ") is not a valid choice for RetVal.
         RetVal needs to be 0 (p-value) or 1 (test statistic) or 2 (test statistic an p-value).")
  }
  
  if (sum(is.na(freq))>0){
    stop("Allele frequency matrix has missing values")
  }
  
  if (sum(is.na(coverage))>0){
    stop("Coverage matrix has missing values")
  }
  
  if(sum(freq<0)>0) {
    stop("Negative allele frequencies are not allowed")
  }
  
  if(sum(coverage<=0)>0) {
    stop("Negative and 0 coverages are not allowed")
  }
  
  if(!is.null(poolSize) && sum(poolSize<=0)>0) {
    stop("Negative and 0 sizes for pool size are not allowed")
  }
  
  if (TA==TRUE && IntGen==FALSE){
    warning("Is it only possible to do Taylor approximation of the variance of the test 
            statistic if intermediate generatons are given. TA option will be ignored")
  }
  
  if (IntGen==FALSE && length(unique(gen))>2){
    stop("IntGen is set to FALSE. Only two generations can be considered.")
  }
  
  if (IntGen==TRUE && length(unique(gen))==2){
    warning("IntGen is set to TRUE, but only two time points are given. They will be used
            as first and last time point and no intermediate generations will be considered.")
    IntGen=FALSE
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
  
  if (length(unique(repl))<2)
    stop("Data for at least 2 replicates are needed")
  
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
  if(npop %% length(gen) != 0 || npop %% length(repl) != 0) 
    stop("The number of populations (", npop, ") has to be a multiple of the length of 'gen'(", length(gen), ") and 'repl' (", length(repl), ").")
  
  ng<-length(unique(gen))
  nreps<-length(unique(repl))
  
  if (length(unique(gen))>1){
    if (order==0) { 
      popInfo <- data.table(pop=1:ncol(freq), gen=c(rep(unique(gen),nreps)), repl=rep(unique(repl),each=ng))
    } else if (order==1) { 
      popInfo <- data.table(pop=1:ncol(freq), gen=c(rep(unique(gen),each=nreps)), repl=rep(unique(repl),ng)) 
    } else {
      stop("The order of the columns in the matrix of allele frequency and coverages can assume only values 0 and 1")
    }
  } else {
    if (order==0) { 
      popInfo <- data.table(pop=1:ncol(freq), gen=c(rep(unique(gen),nreps)), repl=rep(unique(repl),each=2))
    } else if (order==1) { 
      popInfo <- data.table(pop=1:ncol(freq), gen=c(rep(unique(gen),nreps)), repl=rep(unique(repl),2)) 
    } else {
      stop("The order of the columns in the matrix of allele frequency and coverages can assume only values 0 and 1")
    }
  }
  
  
  mask <- rowSums(coverage < mincov)
  stat <- NA
  # no drift
  if(length(unique(popInfo$gen)) == 1) {
    
    if(!missing(Ne))
      warning("Value of 'Ne' will be ignored because no random genetic drift is assumed.")
    
    # individual sequencing
    if(is.null(poolSize)) {
      if (order==0) {
        x1 <- coverage[,seq(1,2*nreps,by=2)]
        x2 <- coverage[,seq(2,2*nreps,by=2)]
        x11 <- freq[,seq(1,2*nreps,by=2)] * x1
        x21 <- freq[,seq(2,2*nreps,by=2)] * x2
      } else {
        x1 <- coverage[,1:nreps]
        x2 <- coverage[,(nreps+1):(2*nreps)]
        x11 <- freq[,1:nreps] * x1
        x21 <- freq[,(nreps+1):(2*nreps)] * x2
      }
      
      n <- x1 + x2
      if(sum(x11==0)>0) {
        x11[x11==0] <- rep(1,sum(x11==0))
        warning('The counts that equal 0 or equal the coverage of the considered 
                locus are changed to 1 or to coverage-1 respectively.')
      }
      if(sum(x11==x1)>0){
        x11[x11==x1] <- x1[x11==x1]-1
        warning('The counts assuming values 0 or equal the coverage of the considered 
                locus are changed to 1 and to coverage-1 respectively.')
      } 
      x12 <- x1 - x11
      x22 <- x2 - x21
      x_1 <- x11 + x21
      x_2 <- x12 + x22
      if (nrow(freq)==1) {stat<-(sum(x11 - x1*x_1/n))^2/sum(x1*x_1*x2*x_2/(n^2*(n-1))) }
      else { stat<-(rowSums(x11 - x1*x_1/n))^2/rowSums(x1*x_1*x2*x_2/(n^2*(n-1)))  }
      # pooled sequencing
      } else if(npop %% length(poolSize) == 0) {
        if (order==0) {
          R1 <- coverage[,seq(1,2*nreps,by=2)]
          R2 <- coverage[,seq(2,2*nreps,by=2)]
          x11R <- freq[,seq(1,2*nreps,by=2)] * R1
          x21R <- freq[,seq(2,2*nreps,by=2)] * R2
          x1 <- poolSize[seq(1,2*nreps, by=2)]
          x2 <- poolSize[seq(2,2*nreps, by=2)]
        } else {
          R1 <- coverage[,1:nreps]
          R2 <- coverage[,(nreps+1):(2*nreps)]
          x11R <- freq[,1:nreps] * R1
          x21R <- freq[,(nreps+1):(2*nreps)] * R2
          x1 <- poolSize[1:nreps]
          x2 <- poolSize[(nreps+1):(2*nreps)]
        }
        n<- R1+R2
        
        if(sum(x11R==0)>0) {
          x11R[x11R==0] <- rep(1,sum(x11R==0))
          warning('The counts assuming values 0 or equal the coverage of the considered 
                  locus are changed to 1 and to coverage-1 respectively.')
        }
        if(sum(x11R==R1)>0){
          x11R[x11R==R1] <- R1[x11R==R1]-1
          warning('The counts assuming values 0 or equal the coverage of the considered 
                  locus are changed to 1 and to coverage-1 respectively.')
        } 
        x12R <- R1 - x11R
        x22R <- R2 - x21R
        
        x_1<-x11R+x21R
        x_2<-x12R+x22R
        
        
        if (nrow(freq)==1) {stat <- sum(x11R - R1*x_1/n)^2/sum((R2/n)^2*(x11R*x12R/R1*(1+(R1-1)/x1)) + (R1/n)^2*(x21R*x22R/R2*(1+(R2-1)/x2))) }
        else { stat <- rowSums(x11R - R1*x_1/n)^2/rowSums((R2/n)^2*(x11R*x12R/R1*(1+(R1-1)/x1)) + (R1/n)^2*(x21R*x22R/R2*(1+(R2-1)/x2))) }
        } else {
          stop("The number of populations (", npop, ") has to be a multiple of the length of 'poolSize' (", length(poolSize), ").")
        }
    
    # with drift #####
      } else { 
        if (length(Ne)!=length(unique(repl))){
          stop("For each replicate a corresponding value of the effective population size 'Ne' has to be given")
        }
        
        # individual sequencing
        if(is.null(poolSize)) {
          ming <- min(gen)
          maxg <- max(gen)
          x1 <- coverage[,popInfo$gen == ming]
          x2 <- coverage[,popInfo$gen == maxg]
          n<-x1+x2
          
          x11 <- matrix(freq[,popInfo$gen == ming], nrow = nrow(freq)) * x1
          
          if(sum(x11==0)>0) {
            x11[x11==0] <- rep(1,sum(x11==0))
            warning('The counts assuming values 0 or equal the coverage of the considered 
                    locus are changed to 1 and to coverage-1 respectively.')
          }
          if(sum(x11==x1)>0){
            x11[x11==x1] <- x1[x11==x1]-1
            warning('The counts assuming values 0 or equal the coverage of the considered 
                    locus are changed to 1 and to coverage-1 respectively.')
          } 
          
          x12 <- x1 - x11
          x21 <- matrix(freq[,popInfo$gen == maxg], nrow = nrow(freq)) * x2
          x22 <- x2 - x21
          x_1<-x11+x21
          x_2<-x12+x22
          if (nrow(freq)==1) {sigd <- x11/x1*(1-x11/x1)*(1-(1-1/(2*Ne))^(maxg-ming))}
          else { sigd <- x11/x1*(1-x11/x1)*matrix(rep(1-(1-1/(2*Ne))^(maxg-ming),nrow(x1)), nrow = nrow(x1), byrow = TRUE) }
          
          if (IntGen == TRUE){
            sigd <- c()
            freq_int <- matrix(freq[,!(popInfo$gen %in% c(ming, maxg))], nrow = nrow(freq))
            coverage_int <- matrix(coverage[,!(popInfo$gen %in% c(ming, maxg))], nrow = nrow(freq))
            
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
            
            x_int2 <- coverage_int - x_int1
             if (order==0){
              if (TA==FALSE){
                if (nrow(freq)==1) { sigd_mat <- (c(x11,x_int1)*(c(x1,coverage_int)-c(x11,x_int1))/c(x1,coverage_int)^2)*((1-(1-1/(2*rep(Ne,each=(ng-1))))^(c(rep(sort(gen)[2]-ming,each=nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],nreps)))))}
                else { sigd_mat <- (cbind(x11,x_int1)*(cbind(x1,coverage_int)-cbind(x11,x_int1))/cbind(x1,coverage_int)^2)*(matrix(rep(1-(1-1/(2*rep(Ne,each=(ng-1))))^(c(rep(sort(gen)[2]-ming,each=nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],nreps))),nrow(x1)),nrow = nrow(x1), byrow = TRUE)) }
              } 
              if (TA==TRUE){
                if (nrow(freq)==1) { sigd_mat <- (c(x11,x_int1)*(c(x1,coverage_int)-c(x11,x_int1))/c(x1,coverage_int)^2)*(c(rep(sort(gen)[2]-ming,each=nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],nreps)))/((2*rep(Ne,each=(ng-1)))) }
                else { sigd_mat <- (cbind(x11,x_int1)*(cbind(x1,coverage_int)-cbind(x11,x_int1))/cbind(x1,coverage_int)^2)*(matrix(rep(c(rep(sort(gen)[2]-ming,each=nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],nreps))/((2*rep(Ne,each=(ng-1)))),nrow(x1)),nrow = nrow(x1), byrow = TRUE)) }
              }
            } else {
              if (TA==FALSE){
                if (nrow(freq)==1) { sigd_mat <- (c(x11,x_int1)*(c(x1,coverage_int)-c(x11,x_int1))/c(x1,coverage_int)^2)*((1-(1-1/(2*rep(Ne,ng-1)))^(c(rep(sort(gen)[2]-ming,nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],each=nreps)))))}
                else { sigd_mat <- (cbind(x11,x_int1)*(cbind(x1,coverage_int)-cbind(x11,x_int1))/cbind(x1,coverage_int)^2)*(matrix(rep(1-(1-1/(2*rep(Ne,ng-1)))^(c(rep(sort(gen)[2]-ming,nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],each=nreps))),nrow(x1)),nrow = nrow(x1), byrow = TRUE)) }
              } 
              if (TA==TRUE){
                if (nrow(freq)==1) { sigd_mat <- (c(x11,x_int1)*(c(x1,coverage_int)-c(x11,x_int1))/c(x1,coverage_int)^2)*(c(rep(sort(gen)[2]-ming,nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],each=nreps)))/((2*rep(Ne,ng-1))) }
                else { sigd_mat <- (cbind(x11,x_int1)*(cbind(x1,coverage_int)-cbind(x11,x_int1))/cbind(x1,coverage_int)^2)*(matrix(rep(c(rep(sort(gen)[2]-ming,nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],each=nreps)),nrow(x1))/((2*rep(Ne,ng-1))),nrow = nrow(x1), byrow = TRUE)) }
              }
            }
            

            for (r in unique(repl)){
              if (is.vector(sigd_mat)==TRUE) {sigd <- c(sigd, sum(sigd_mat[c(popInfo$repl[popInfo$gen==ming]==r,popInfo$repl[!(popInfo$gen %in% c(ming, maxg))]==r)]))}
              else {sigd <- cbind(sigd,rowSums(sigd_mat[,c(popInfo$repl[popInfo$gen==ming]==r,popInfo$repl[!(popInfo$gen %in% c(ming, maxg))]==r)]))}
            }
            }

          Ep2 <- (x11/x1 + x21/x2)/2
          if(MeanStart==FALSE)
            Ep2 <- x11/x1
          
          if (nrow(freq)==1) {stat <- sum(x11 - x1*x_1/n)^2/sum((x2/n)^2*(x11*x12/x1) + (x1/n)^2*(x2*Ep2*(1-Ep2)+x2*(x2-1)*(sigd)))}
          else { stat <- rowSums(x11 - x1*x_1/n)^2/rowSums((x2/n)^2*(x11*x12/x1) + (x1/n)^2*(x2*Ep2*(1-Ep2)+x2*(x2-1)*(sigd)))  }

          } else if(npop %% length(poolSize) == 0) {
            ming <- min(gen)
            maxg <- max(gen)
            R1 <- matrix(coverage[,popInfo$gen == ming], nrow = nrow(coverage))
            R2 <- matrix(coverage[,popInfo$gen == maxg], nrow = nrow(coverage))
            n <- R1 + R2
            
            x11R <- matrix(freq[,popInfo$gen == ming], nrow = nrow(freq)) * R1
            
            if(sum(x11R==0)>0) {
              x11R[x11R==0] <- rep(1,sum(x11R==0))
              warning('The counts assuming values 0 or equal the coverage of the considered 
                      locus are changed to 1 and to coverage-1 respectively.')
            }
            if(sum(x11R==R1)>0){
              x11R[x11R==R1] <- R1[x11R==R1]-1
              warning('The counts assuming values 0 or equal the coverage of the considered 
                      locus are changed to 1 and to coverage-1 respectively.')
            } 
            
            x12R <- R1 - x11R
            
            x21R <- matrix(freq[,popInfo$gen == maxg], nrow = nrow(freq)) * R2
            x22R <- R2 - x21R
            
            x1 <- poolSize[popInfo$gen == ming]
            x2 <- poolSize[popInfo$gen == maxg]
            
            x_1 <- x11R+x21R
            x_2 <- x12R+x22R
            sigd <- x11R/R1*(1-x11R/R1)*(1-(1-1/(2*Ne))^(maxg-ming))
            
            if (IntGen == TRUE){
              sigd <- c()
              freq_int <- matrix(freq[,!(popInfo$gen %in% c(ming, maxg))], nrow = nrow(freq))
              coverage_int <- matrix(coverage[,!(popInfo$gen %in% c(ming, maxg))], nrow = nrow(freq))
              
              
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
              
              x_int2 <- coverage_int - x_int1
              pool_int <- coverage_int
              
              
              if (order==0){
                if (TA==FALSE){
                  if (nrow(freq)==1) {sigd_mat <- (c(x11R,x_int1)*(c(R1,pool_int)-c(x11R,x_int1))/c(R1,pool_int)^2)*((1-(1-1/(2*rep(Ne,each=(ng-1))))^(c(rep(sort(gen)[2]-ming,each=nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],nreps)))))}
                  else { sigd_mat <- (cbind(x11R,x_int1)*(cbind(R1,pool_int)-cbind(x11R,x_int1))/cbind(R1,pool_int)^2)*(matrix(rep(1-(1-1/(2*rep(Ne,each=(ng-1))))^(c(rep(sort(gen)[2]-ming,each=nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],nreps))),nrow(R1)),nrow = nrow(R1), byrow = TRUE)) }
                } 
                if (TA==TRUE) {
                  if (nrow(freq)==1) { sigd_mat <- (c(x11R,x_int1)*(c(R1,pool_int)-c(x11R,x_int1))/c(R1,pool_int)^2)*(c(rep(sort(gen)[2]-ming,each=nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],nreps)))/((2*rep(Ne,each=(ng-1))))}
                  else { sigd_mat <- (cbind(x11R,x_int1)*(cbind(R1,pool_int)-cbind(x11R,x_int1))/cbind(R1,pool_int)^2)*(matrix(rep(c(rep(sort(gen)[2]-ming,each=nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],nreps)),nrow(R1))/((2*rep(Ne,each=(ng-1)))),nrow = nrow(R1), byrow = TRUE)) }
                }
              } else {
                if (TA==FALSE){
                  if (nrow(freq)==1) {sigd_mat <- (c(x11R,x_int1)*(c(R1,pool_int)-c(x11R,x_int1))/c(R1,pool_int)^2)*((1-(1-1/(2*rep(Ne,ng-1)))^(c(rep(sort(gen)[2]-ming,nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],each=nreps)))))}
                  else { sigd_mat <- (cbind(x11R,x_int1)*(cbind(R1,pool_int)-cbind(x11R,x_int1))/cbind(R1,pool_int)^2)*(matrix(rep(1-(1-1/(2*rep(Ne,ng-1)))^(c(rep(sort(gen)[2]-ming,nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],each=nreps))),nrow(R1)),nrow = nrow(R1), byrow = TRUE)) }
                } 
                if (TA==TRUE) {
                  if (nrow(freq)==1) { sigd_mat <- (c(x11R,x_int1)*(c(R1,pool_int)-c(x11R,x_int1))/c(R1,pool_int)^2)*(c(rep(sort(gen)[2]-ming,nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],each=nreps)))/((2*rep(Ne,ng-1)))}
                  else { sigd_mat <- (cbind(x11R,x_int1)*(cbind(R1,pool_int)-cbind(x11R,x_int1))/cbind(R1,pool_int)^2)*(matrix(rep(c(rep(sort(gen)[2]-ming,nreps),rep(sort(gen)[-c(1,2)]-sort(gen)[-c(1,ng)],each=nreps)),nrow(R1))/((2*rep(Ne,ng-1))),nrow = nrow(R1), byrow = TRUE)) }
                }
              }
              
              for (r in unique(repl)){
                if (is.vector(sigd_mat)==TRUE) {sigd <- c(sigd, sum(sigd_mat[c(popInfo$repl[popInfo$gen==ming]==r,popInfo$repl[!(popInfo$gen %in% c(ming, maxg))]==r)]))}
                else {sigd <- cbind(sigd,rowSums(sigd_mat[,c(popInfo$repl[popInfo$gen==ming]==r,popInfo$repl[!(popInfo$gen %in% c(ming, maxg))]==r)]))}
              }
              }
            
            
            Ep2 <- (x11R/R1 + x21R/R2)/2 
            if(MeanStart==FALSE){Ep2 <- x11R/R1}
            
            if (nrow(freq)==1) {stat <- sum(x11R - R1*x_1/n)^2/sum((R2/n)^2*(x11R*x12R/R1*(1+(R1-1)/x1)) + (R1/n)^2*(R2*Ep2-R2*Ep2^2+(R2*(R2-1)/x2)*(Ep2*(1-Ep2)+(x2-1)*(sigd))))  }
            else { stat <- rowSums(x11R - R1*x_1/n)^2/rowSums((R2/n)^2*(x11R*x12R/R1*(1+(R1-1)/x1)) + (R1/n)^2*(R2*Ep2-R2*Ep2^2+(R2*(R2-1)/matrix(rep(x2,nrow(R2)),nrow=nrow(R2)))*(Ep2*(1-Ep2)+(matrix(rep(x2-1,nrow(Ep2)),nrow = nrow(Ep2)))*(sigd)))) }
            
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
