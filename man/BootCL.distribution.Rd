\name{BootCL.distribution}
\alias{BootCL.distribution}

\title{Find the bootstrap distribution}
\description{
 BootCL.distribution finds the bootstrap distribution for the
     measure of significance of  chromosomal clustering in
     differentially expressed genes defined from a microarray
     experiment.
}
\usage{
 BootCL.distribution(chip, sampling.count, total.sampling.count,windowsize)
}

\arguments{
  \item{chip}{name of chip}
  \item{sampling.count}{the number of genes in the cluster}
  \item{total.sampling.count}{ the number of bootstrap sampling}
  \item{windowsize}{window size}
}

\value{
  \item{dist}{Distribution from Bootstrap samples}
  \item{windowsize}{window size}
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
}
\keyword{distribution}% at least one, from doc/KEYWORDS
