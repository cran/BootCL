\name{BootCL.plot}
\alias{BootCL.plot}

\title{Display the bootstrap distribution and p-value }
\description{
     BootCL.plot presents the bootstrap distribution and p-value for
     the test  of chromosomal clustering in differentially expressed
     genes defined  from a microarray experiment
}
\usage{
BootCL.plot(BootP,xrange=NULL,freq.bootCL=NULL)
}

\arguments{
  \item{BootP}{The result from BootCL.distribution}
  \item{xrange}{the range of x-axis}
  \item{freq.bootCL}{logical; if 'TRUE', the histogram graphic is a representation
          of frequencies, the 'counts' component of the result; if
          'FALSE', probability densities, component 'density', are
          plotted }
}

\value{
  \item{BootP}{The result from BootCL.distribution}
}
%\references{ ~put references to the literature/web site here ~ }
\examples{
     library(BootCL)
     chip.name<-"HG.U133A"
     data(affy.hgu133a)
     sample.cluster<-affy.hgu133a[sample(1:1000,100),1]
     sample.ID<-find.ID(chip.name,sample.cluster,affyID.flag=TRUE) 
     Boot.Stat<-BootCL.Statistic(chip.name,sample.ID)
     Bdist<-BootCL.distribution(chip.name,Boot.Stat$sampling.count,total.sampling.count=10000)
     BootP<-BootCL.Pvalue(Boot.Stat,Bdist)
     BootCL.plot(BootP)
}
\keyword{distribution}% at least one, from doc/KEYWORDS
