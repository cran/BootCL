\name{find.ID}
\alias{find.ID}
\title{ Match affyID to Accession IDs}
\description{
 Match affyID to Accession IDs that have chromosome information
}
\usage{
find.ID(chip, sample,affyID.flag)
}
\arguments{
  \item{chip}{ name of chip}
  \item{sample}{ affyIDs or Accession IDs in the cluster}
  \item{affyID.flag} {if affyIDs are used for sample, affyID.flag=TRUE. Else if Accession IDs are used , affyID.flag=FALSE}
}
\value{
  \item{Access.list} {Accession IDs}
  \item{Affy.list}{corresponding affyIDs}
}
%\references{ ~put references to the literature/web site here ~ }

\examples{
     library(BootCL)
     chip.name<-"HG.U133A"
     data(affy.hgu133a)
     sample.cluster<-affy.hgu133a[sample(1:1000,100),1]
     sample.ID<-find.ID(chip.name,sample.cluster) 
}
\keyword{distribution }
