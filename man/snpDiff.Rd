\name{snpDiff}
\alias{snpDiff}
\title{
main function to test point differences between the two bam files
}
\description{
test base count (A,T,G,C) difference at a given position between the two bam files. The base count matrix is tested
by fisher exact test and euclidean distance. If nothing difference, NULL will be returned.
}
\usage{
snpDiff(bam1, bam2, refFsa, regChr, regStart, regEnd, minBaseQuality = 13, minMapQuality = 0, nCores = 1, pValueCutOff = 0.05, baseDistCutOff = 0.1,verbose = TRUE)
}
\arguments{
  \item{bam1}{
the first bam file to be compared
}
  \item{bam2}{
the second bam file to be compared
}
  \item{refFsa}{
the reference fasta file used for bam1 and bam2 alignments
}
  \item{regChr}{
chromosome name of the region of interest, it should match the chromosome name in reference name 
}
  \item{regStart}{
the start position (1-based) of the region of interest
}
  \item{regEnd}{
the end position (1-based) of the region of interest
}
  \item{minBaseQuality}{
the minimum base quality to be used for base count
}
  \item{minMapQuality}{
the minimum read mapping quality to be used for base count
}
  \item{nCores}{
number of thread used for calculate in parallel
}
  \item{pValueCutOff}{
p.value cutoff from fisher.test to display output. If there is no difference between two compared positions (p.value = 1 and d.value = 0), NULL will be returned even setting pValueCutOff = 1.
}
  \item{baseDistCutOff}{
euclidean distance cutoff from dist(,method='euclidean') to display output. If there is no difference between two compared positions (p.value = 1 and d.value = 0), NULL will be returned even setting baseDistCutOff = 0.
}
  \item{verbose}{
print progress on screen, default = TRUE.
  }
}

\value{\code{snpDiff}: returns a data.frame with difference information: chromosome, position, reference base, base count (A,C,G,T,N) for two bam files, p.value (fisher exact test of these base counts) and d.value (euclidean distance of these base counts)

}

\author{
 Xiaobin Xing, <email:xiaobinxing0316@gmail.com>
 }
 
 \references{
 Morgan M, Pages H, Obenchain V and Hayden N. Rsamtools: Binary alignment (BAM), FASTA, variant call (BCF), and tabix file import.
 }
 
\examples{

bam1 <- system.file(package='SICtools','extdata','example1.bam')
bam2 <- system.file(package='SICtools','extdata','example2.bam')
refFsa <- system.file(package='SICtools','extdata','example.ref.fasta')

snpDiff(bam1,bam2,refFsa,'chr04',962501,1026983,pValueCutOff=1,baseDistCutOff=0)

}

