\name{BootCL.Statistic}
\alias{BootCL.Statistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate a measure of chromosomeal spatial bias}
\description{
     Calculate a measure of chromosomeal spatial bias,
     total number of counts that how many two consecutive genes are located
     in the same chromosome and within a predefined windowsize.
}
\usage{
BootCL.Statistic(chip, ID.data, windowsize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chip}{ name of chip}
  \item{ID.data}{ IDs from find.ID function. }
  \item{windowsize}{ window size}
}

\value{
  \item{Diff.count}{ statistic to test statistical significance of  chromosomal clustering }
  \item{windowsize}{window size}
  \item{sampling.count}{the number of genes in the cluster}
  \item{conseq.state}{indicator whether two consecutive genes are located within a predefined windowsize}
  \item{Affy.ID}{Affymetrix ID}
  \item{Access.ID}{Accession ID}
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

}
\keyword{distribution}
