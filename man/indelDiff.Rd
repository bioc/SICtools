\name{indelDiff}
\alias{indelDiff}
\title{
main function to call indel difference between the bam files
}
\description{
test indel-read count differences at a given indel position between the two bam files. The indel position are obtained
by samtools+bcftools first, and count the number of reads that span no less than 3bp of the indel boundary. The
read-count matrix at a given indel region from the two bam files are tested by fisher exact test and euclidean distance.
If nothing difference, NULL will be returned.
}
\usage{
indelDiff(bam1, bam2, refFsa, regChr, regStart, regEnd, minBaseQuality = 13, minMapQuality = 0, nCores = 1, pValueCutOff= 0.05,gtDistCutOff = 0.1,verbose = TRUE)
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
the minimum base quality to be used for indel-read count
}
  \item{minMapQuality}{
the minimum read mapping quality to be used for indel-read count
}
  \item{nCores}{
number of thread used for calculate in parallel
}
  \item{pValueCutOff}{
p.value cutoff from fisher.test to display output. If there is no difference between two compared positions (p.value = 1 and d.value = 0), NULL will be returned even setting pValueCutOff = 1.
}
  \item{gtDistCutOff}{
euclidean distance cutoff from dist(,method='euclidean') to display output. If there is no difference between two compared positions (p.value = 1 and d.value = 0), NULL will be returned even setting gtDistCutOff = 0.
}
  \item{verbose}{
print progress on screen, default = TRUE.
  }

 }

 \value{\code{indelDiff}: returns a data.frame with difference information: chromosome, position, reference genenotype, two alt genotypes, and their indel-read count for two bam files, p.value (fisher exact test of these read counts) and d.value (euclidean distance of these read counts)
}

 \author{
 Xiaobin Xing, <email:xiaobinxing0316@gmail.com>
 }
 
 \references{
Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. 
and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. 
Bioinformatics, 25, 2078-9. [PMID: 19505943] }
 
 \examples{
bam1 <- system.file(package='SICtools','extdata','example1.bam')
bam2 <- system.file(package='SICtools','extdata','example2.bam')
refFsa <- system.file(package='SICtools','extdata','example.ref.fasta')

indelDiff(bam1,bam2,refFsa,'chr07',828514,828914,pValueCutOff=1,gtDistCutOff=0)

}

