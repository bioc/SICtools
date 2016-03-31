# TODO: test snpDiff() using RUnit

## for normal test
test_snpDiff <- function(){
	
	bam1 <- system.file(package='SICtools','extdata','example1.bam')
	bam2 <- system.file(package='SICtools','extdata','example2.bam')
	refFsa <- system.file(package='SICtools','extdata','example.ref.fasta')
	
	## check
	snpTest <- snpDiff(bam1,bam2,refFsa,'chr04',962501,1026983,pValueCutOff=1,baseDistCutOff=0,nCores=2)
	checkEquals(nrow(snpTest),16)
	checkEquals(subset(snpTest,p.value < 1e-10)$pos,c(962801,1026683))
	
	## check NULL
	checkTrue(is.null(snpDiff(bam1,bam2,refFsa,'chr04',962501,1026983,pValueCutOff=0,baseDistCutOff=0,nCores=2)))
	
}

