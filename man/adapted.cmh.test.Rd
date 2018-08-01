\name{adapted.cmh.test}
\alias{adapted.cmh.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adapted Cochran-Mantel-Haenszel (CMH) test to drift and pool sequencing variance. 
}
\description{
Cochran-Mantel-Haenszel test can be used to detect the presence of selection in time series data with replicates populations. The function adapted.cmh.test performs a modified version of the classical CMH test (Agresti 2002). The variance of the test statistic is modified in order to account for the different component of the variance in genetic time series data. Indeed, the variance term now includes the drift variance and the variance introduced by pool sequencing. The output of the function is a vector which contains one p-value (and/or one value of the test statistic) for each SNP provided in the input data.
}
\usage{
adapted.cmh.test(freq, coverage, Ne, gen, repl, poolSize = NULL, 
mincov = 1, MeanStart = TRUE, IntGen = FALSE, TA = FALSE, order = 0, 
RetVal = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{freq}{
numeric matrix containing the allele frequency information. Each row corresponds to a SNP and each column corresponds to a specific replicate in a specific time point.
}
  \item{coverage}{
numeric matrix containing the coverage information for each SNP, for each replicate population at each time point. Same structure of freq.
}
  \item{Ne}{
numeric vector specifying the known or estimated value for the effective population size for each replicate population. The order of the replicates must respect the structure of freq and coverage.
}
  \item{gen}{
numeric vector containing the time points in which the allele freqeuncy data are recorded. The order of the generations must respect the structure of freq and coverage.
}
  \item{repl}{
numeric vector containing the replicates for which the allele freqeuncy data are recorded. The order of the replicates must respect the structure of freq and coverage.
}

  \item{poolSize}{
numeric vector containing the value of the sample size used for pool sequencing in each time point and for each replicate. The length and the structure of the vector must respect the one in freq and coverage. If no pool sequencing step occurs then poolSize = NULL.
}
  \item{mincov}{
numeric specyfing the allowed minimum value for the coverage. When this requirement is not fulfilled the test statistic for that specific SNP is not computed and NA are returned. 
}
  \item{MeanStart}{
logical. If MeanStart = FALSE the read counts of the first time point are used to compute the expected allele frequency of the last time point in the variance estimator. If MeanStart = TRUE also the intermediate and the last time points are used for this purpose and the estimator corresponds to the average of the counts of all the time points.
}
  \item{IntGen}{
logical. IntGen = TRUE when the allele frequency and coverage matrix with information of the intermediate time points is provided, IntGen = FALSE otherwise.
}
  \item{TA}{
logical. TA = TRUE if Taylor approximation for the variance computation has to be used. TA = FALSE otherwise.
}
  \item{order}{
factor. Suppose to have N generations and K replicates. If order = 0 then the columns in freq are assumed to follow the structure: generation0.replicate1, generation1.replicate1, generation2.replicate1, ..., generationN.replicate1, generation0.replicate2, ..., generationN.replicate2, ... ..., generation0.replicateK, ..., generationN.replicateK. If order = 1 then the columns in freq are assumed to follow the structure: generation0.replicate1, generation0.replicate2, generation0.replicate3, ..., generation0.replicateK, generation1.replicate1, ..., generation1.replicateK, ... ..., generationN.replicate1, ..., generationN.replicateK. 
}
  \item{RetVal}{
factor. If RetVal = 0 a vector containing one p-value per SNP is returned. If RetVal = 1 a vector containing one value of the test statistic per SNP is returned. If RetVal = 2 a matrix containing the test statistic values in the first column and the p-values in the second column is returned.
}
}

\value{
adapted.cmh.test returns either a vector (RetVal = 0 or RetVal = 1) or a matrix (RetVal = 2). If RetVal = 0 or RetVal = 1 then the length of the vector equals the number of rows of the freq matrix of input and it contains the p-values (RetVal = 0) or the test statistic values (RetVal = 1) for each SNP. If RetVal = 2 adapted.cmh.test returns a matrix with 2 columns and number of rows equal to the number of rows of the freq matrix of input. Here the first column contains the values of the test statistic and the second column the p-values for each SNP.
}
\references{
Alan Agresti. Categorical Data Analysis. Wiley Interscience, New York, 2002. ISBN 0-471-36093-7 978-0-471-36093-3.
}
\author{
Marta Pelizzola, Kerstin Gaertner
}

\seealso{
\code{\link{adapted.chisq.test}}.
}
\examples{
#The following examples use simulated data. In order to use 
#use the function with real data stored in the common sync 
#file format we suggest to use the R package poolSeq 
#(https://github.com/ThomasTaus/poolSeq) in
#order to store them in the proper format required from this 
#function

##########################################
############Classical CMH test############
##########################################

#simulate allele frequency for 3 replicates and coverage matrix
rep <- c(1,2,3)
p <- matrix(c(rep(runif(10000, min=0, max=1),2),rep(runif(10000, 
            min=0, max=1),2),rep(runif(10000, min=0, max=1),2)),
            ncol=length(rep), byrow = FALSE)
covMat <- cbind(rpois(nrow(p), lambda=200), rpois(nrow(p), lambda=80), 
                rpois(nrow(p), lambda=200), rpois(nrow(p), lambda=80), 
                rpois(nrow(p), lambda=200), rpois(nrow(p), lambda=80))

#add sampling noise to the allele frequency matrix
afMat <- rbinom(n=length(covMat), size=covMat, prob=p) / covMat

#compute the appropriate CMH test for this situation: no drift 
#(i.e. only one time point) and no pool sequencing (i.e.
#poolSize = NULL)
pval <- adapted.cmh.test(freq=afMat, coverage=covMat, gen=0, 
                                                repl = rep)

###################################################################
################CMH test adapted to pool sequencing################
###################################################################

#simulate allele frequency and coverage matrix
rep <- c(1,2,3)
p <- matrix(c(rep(runif(10000, min=0, max=1),2),rep(runif(10000, 
            min=0, max=1),2),rep(runif(10000, min=0, max=1),2)),
            ncol=length(rep), byrow = FALSE)
covMat <- cbind(rpois(length(p), lambda=80), rpois(length(p), lambda=80), 
                rpois(length(p), lambda=80), rpois(length(p), lambda=80), 
                rpois(length(p), lambda=80), rpois(length(p), lambda=80))

#add pool sequencing noise
ps <- rep(1000, 2*length(rep))
afMat <- rbinom(n=length(covMat), size=matrix(rep(ps,nrow(covMat)), 
          ncol=2*length(rep)), prob=p) / matrix(rep(ps,nrow(covMat)), 
          ncol=2*length(rep))

#add sampling noise to the allele frequency matrix
afMat <- rbinom(n=length(covMat), size=covMat, prob=afMat) / covMat

#compute the appropriate CMH test for this situation: no drift (i.e. 
#only one time point) and pool sequencing (i.e.
#poolSize = 1000)
pval <- adapted.cmh.test(freq=afMat, coverage=covMat, gen=0, 
                                    repl = rep, poolSize=ps)


\dontrun{

## NOTE: need poolSeq (https://github.com/ThomasTaus/poolSeq) to 
#simulate drift in the following two examples:

###################################################################
#####################CMH test adapted to drift#####################
###################################################################

#simulate the starting allele frequencies 
rep <- c(1,2,3)
p <- matrix(c(round(runif(10000, min=0, max=1)*2*Ne)/(2*Ne),
              round(runif(10000, min=0, max=1)*2*Ne)/(2*Ne),
              round(runif(10000, min=0, max=1)*2*Ne)/(2*Ne)), 
              ncol=length(rep))

#simulate random drift with generations 0 and 60
Ne <- 300
tp <- c(0, 60)
p_ev1 <- wf.traj(p0=p[,1], Ne=Ne, t=tp)
p_ev2 <- wf.traj(p0=p[,2], Ne=Ne, t=tp)
p_ev3 <- wf.traj(p0=p[,3], Ne=Ne, t=tp)
afMat <- cbind(p_ev1, p_ev2, p_ev3)

#simulate the coverage matrix
covMat <- cbind(rpois(length(p), lambda=80), rpois(length(p), lambda=80), 
                rpois(length(p), lambda=80), rpois(length(p), lambda=80), 
                rpois(length(p), lambda=80), rpois(length(p), lambda=80))

#add sampling noise to the allele frequency matrix
afMat <- rbinom(n=length(covMat), size=covMat, prob=afMat) / covMat

#compute the appropriate CMH test for this situation: drift (i.e. 
#evolving population for 60 generations) and no pool sequencing 
#(i.e. poolSize = NULL)
pval <- adapted.cmh.test(freq=afMat, coverage=covMat, Ne=rep(Ne,
                                length(rep)), gen=tp, repl = rep)


###########################################################################
################CMH test adapted to drift & pool sequencing################
###########################################################################


#simulate the starting allele frequencies 
rep <- c(1,2,3)
p <- matrix(c(round(runif(10000, min=0, max=1)*2*Ne)/(2*Ne),
            round(runif(10000, min=0, max=1)*2*Ne)/(2*Ne),
            round(runif(10000, min=0, max=1)*2*Ne)/(2*Ne)), 
            ncol=length(rep))

#simulate random drift with generations 0 and 60
Ne <- 300
tp <- c(0, 60)
p_ev1 <- wf.traj(p0=p[,1], Ne=Ne, t=tp)
p_ev2 <- wf.traj(p0=p[,2], Ne=Ne, t=tp)
p_ev3 <- wf.traj(p0=p[,3], Ne=Ne, t=tp)
afMat <- cbind(p_ev1, p_ev2, p_ev3)

#add pool sequencing noise
ps <- rep(1000, 2*length(rep))
afMat <- rbinom(n=length(covMat), size=matrix(rep(ps,nrow(covMat)), 
          ncol=2*length(rep)), prob=p) / matrix(rep(ps,nrow(covMat)), 
          ncol=2*length(rep))

#add sampling noise to the allele frequency matrix
afMat <- rbinom(n=length(covMat), size=covMat, prob=afMat) / covMat

#compute the appropriate CMH test for this situation: drift 
#(i.e. evolving population for 60 generations) and pool 
#sequencing (i.e. poolSize = 1000)
pval <- adapted.cmh.test(freq=afMat, coverage=covMat, Ne=rep(Ne,
                    length(rep)), gen=tp, repl=rep, poolSize=ps)
    }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{selection}% use one of  RShowDoc("KEYWORDS")
%\keyword{time series data}% __ONLY ONE__ keyword per line
%\keyword{replicate populations}