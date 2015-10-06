\name{SICtools-package}
\alias{SICtools-package}
\alias{SICtools}
\docType{package}
\title{
tools for SNV/Indel Comparison between two bam files with near relationship
}
\description{
             This package is to find SNV/Indel differences between two bam files with near relationship
             in a way of pairwise comparison thourgh each base position across the genome region of interest. 
             The difference is inferred by fisher test and euclidean distance, the input of which is the base
             count (A,T,G,C) in a given position and read counts for indels that span no less than 2bp on both 
             sides of indel region called from samtools+bcftools
}
\details{
\tabular{ll}{
Package: \tab SICtools\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-07-24\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab Yes\cr
}

}
\author{
Xiaobin Xing

Maintainer: Xiaobin Xing <xxb0316@gmail.com>
}

\keyword{ SICtools }
