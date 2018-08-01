\name{adapted.chisq.test}
\alias{adapted.chisq.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adapted Chi-squared test to drift and pool sequencing variance. 
}
\description{
Chi-squared test can be used to detect the presence of selection in time series data. The function adapted.chisq.test performs a modified version of the classical Chi-squared test (Pearson 1900). The variance of the test statistic is modified in order to account for the different component of the variance in genetic time series data. Indeed, the variance term now includes the drift variance and the variance introduced by pool sequencing. The output of the function is a vector which contains one p-value (and possibly one value of the test statistic) for each SNP provided in the input data.
}
\usage{
adapted.chisq.test(freq, coverage, Ne, gen, poolSize = NULL, 
mincov = 1, MeanStart = TRUE, IntGen = FALSE, TA = FALSE, RetVal = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{freq}{
numeric matrix containing the allele frequency information. Each row corresponds to a SNP and each column corresponds to a time point.
}
  \item{coverage}{
numeric matrix containing the coverage information for each SNP at each time point. Same structure of freq.
}
  \item{Ne}{
integer specifying the known or estimated value for the effective population size (numeric).
}
  \item{gen}{
numeric vector containing the time points in which the allele freqeuncy data are recorded. The order of the generations must respect the structure of freq and coverage.
}
  \item{poolSize}{
numeric vector (same length and same order of gen) containing the value of the sample size used for pool sequencing in each time point. If no pool sequencing step occurs then poolSize = NULL.
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
  \item{RetVal}{
factor. If RetVal = 0 a vector containing one p-value per SNP is returned. If RetVal = 1 a vector containing one value of the test statistic per SNP is returned. If RetVal = 2 a matrix containing the test statistic values in the first column and the p-values in the second column is returned.
}
}

\value{
adapted.chisq.test returns either a vector (RetVal = 0 or RetVal = 1) or a matrix (RetVal = 2). If RetVal = 0 or RetVal = 1 then the length of the vector equals the number of rows of the freq matrix of input and it contains the p-values (RetVal = 0) or the test statistic values (RetVal = 1) for each SNP. If RetVal = 2 adapted.chisq.test returns a matrix with 2 columns and number of rows equal to the number of rows of the freq matrix of input. Here the first column contains the values of the test statistic and the second column the p-values for each SNP.
}
\references{
Pearson, Karl. (1900) On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling. Philosophical Magazine. Series 5. 50 (302) 157-175.
}
\author{
Marta Pelizzola, Kerstin Gaertner
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{adapted.cmh.test}}.
}
\examples{
#The following examples use simulated data. In order to use 
#use the function with real data stored in the common sync 
#file format we suggest to use the R package poolSeq 
#(https://github.com/ThomasTaus/poolSeq) in
#order to store them in the proper format required from this 
#function

############################################################
#################Classical Chi-squared test#################
############################################################

#simulate allele frequency and coverage matrix
p <- runif(10000, min=0, max=1)
covMat <- cbind(rpois(length(p), lambda=80), rpois(length(p), 
                                                  lambda=80))

#add sampling noise to the allele frequency matrix
afMat <- rbinom(n=length(covMat), size=covMat, prob=p) / covMat

#compute the appropriate chi-squared test for this situation: no drift 
#(i.e. only one time point) and no pool sequencing (i.e. 
#poolSize = NULL)
pval <- adapted.chisq.test(freq=afMat, coverage=covMat, gen=0)

#############################################################################
#################Chi-squared test adapted to pool sequencing#################
#############################################################################

#simulate allele frequency and coverage matrix
p <- runif(10000, min=0, max=1)
covMat <- cbind(rpois(length(p), lambda=80), rpois(length(p), 
                                                  lambda=80))

#add pool sequencing noise
ps <- rep(1000, 2)
afMat <- rbinom(n=length(covMat), size=matrix(rep(ps,nrow(covMat)), 
            ncol=2), prob=p) / matrix(rep(ps,nrow(covMat)), ncol=2)

#add sampling noise to the allele frequency matrix
afMat <- rbinom(n=length(covMat), size=covMat, prob=afMat) / covMat

#compute the appropriate chi-squared test for this situation: no drift 
#(i.e. only one time point) and pool sequencing (i.e. poolSize = 1000)
pval <- adapted.chisq.test(freq=afMat, coverage=covMat, gen=0, 
                                                  poolSize=ps)

\dontrun{

## NOTE: need poolSeq (https://github.com/ThomasTaus/poolSeq) to 
#simulate drift in the following two examples:

#############################################################################
######################Chi-squared test adapted to drift######################
#############################################################################

#simulate the starting allele frequencies 
p <- round(runif(10000, min=0, max=1))

#simulate random drift with generations 0 and 60
Ne <- 300
tp <- c(0, 60)
afMat <- wf.traj(p0=(p*2*Ne)/(2*Ne), Ne=Ne, t=tp)

#simulate the coverage matrix
covMat <- matrix(rpois(length(afMat), lambda=80), ncol=ncol(afMat))

#add sampling noise to the allele frequency matrix
afMat <- rbinom(n=length(covMat), size=covMat, prob=afMat) / covMat

#compute the appropriate chi-squared test for this situation: drift 
#(i.e. evolving population for 60 generations) and no pool 
#sequencing (i.e. poolSize = NULL)
pval <- adapted.chisq.test(freq=afMat, coverage=covMat, Ne=Ne, gen=tp)


###################################################################################
################Chi-squared test adapted to drift & pool sequencing################
###################################################################################


#simulate the starting allele frequencies 
p <- round(runif(10000, min=0, max=1))

#simulate random drift with generations 0 and 60
Ne <- 300
tp <- c(0, 60)
afMat <- wf.traj(p0=(p*2*Ne)/(2*Ne), Ne=Ne, t=tp)

#simulate the coverage matrix
covMat <- matrix(rpois(length(afMat), lambda=80), ncol=ncol(afMat))

#add pool sequencing noise
ps <- rep(1000, 2)
afMat <- rbinom(n=length(covMat), size=matrix(rep(ps,nrow(covMat)), 
            ncol=2), prob=p) / matrix(rep(ps,nrow(covMat)), ncol=2)

#add sampling noise to the allele frequency matrix
afMat <- rbinom(n=length(covMat), size=covMat, prob=afMat) / covMat

#compute the appropriate chi-squared test for this situation: drift 
#(i.e. evolving population for 60 generations) and pool 
#sequencing (i.e. poolSize = 1000)
pval <- adapted.chisq.test(freq=afMat, coverage=covMat, Ne=Ne, gen=tp, 
                                                          poolSize=ps)
  }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{selection}% use one of  RShowDoc("KEYWORDS")
%\keyword{time series data}% __ONLY ONE__ keyword per line
