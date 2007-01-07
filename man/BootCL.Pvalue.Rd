\name{BootCL.Pvalue}
\alias{BootCL.Pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the P-value}
\description{
     Calculate the P-value to test statistical significance of 
     chromosomal clustering in differentially expressed genes defined 
     from a microarray experiment.
}
\usage{
BootCL.Pvalue(Bstat,distribution)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Bstat}{ object from BootCL.Statistic function}
  \item{distribution}{  Bootstrap distribution from BootCL.distribution function}
}

\value{
  \item{Diff.count}{ statistic to test statistical significance of  chromosomal clustering }
  \item{Pvalue}{ p-value}
  \item{dist}{Bootstrap distribution}
  \item{windowsize}{window size}
  \item{sampling.count}{the number of genes in the cluster}
  \item{Affy.ID}{Affymetrix ID}
  \item{Access.ID}{Accesssion ID}
  ...
}
%\references{ ~put references to the literature/web site here ~ }
\examples{
     library(BootCL)
     chip.name<-"HG.U133A"
     data(affy.hgu133a)
     sample.cluster<-affy.hgu133a[sample(1:1000,100),1]
     sample.ID<-find.ID(chip.name,sample.cluster) 
     Boot.Stat<-BootCL.Statistic(chip.name,sample.ID)
     Bdist<-BootCL.distribution(chip.name,Boot.Stat$sampling.count,total.sampling.count=10000)
     BootP<-BootCL.Pvalue(Boot.Stat,Bdist)
     BootCL.plot(BootP)
}
\keyword{distribution}
