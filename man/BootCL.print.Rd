\name{BootCL.print}
\alias{BootCL.print}

\title{Display gene lists that are within windowsize consecutively  in the same chromosome}
\description{
     BootCL.print presents gene lists that are within consecutive windowsize in the same chromosome.
}
\usage{
BootCL.print(BootP,cutoff=3,affyID.flag=TRUE)
}

\arguments{
  \item{BootP}{The result from BootCL.distribution}
  \item{cutoff}{the size of gene lists}
  \item{affyID.flag} {if affyIDs are used for sample, affyID.flag=TRUE. Else if Accession IDs are used , affyID.flag=FALSE}
}

\value{
  \item{print.list}{all the gene lists}
  \item{select.list}{selected list that the size of gene lists are larger than cutoff}
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
     BootCL.print(BootP)
}
\keyword{distribution}% at least one, from doc/KEYWORDS
