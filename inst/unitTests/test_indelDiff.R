# TODO: test indelDiff()
# 
test_indelDiff <- function(){
	
	bam1 <- system.file(package='SICtools','extdata','example1.bam')
	bam2 <- system.file(package='SICtools','extdata','example2.bam')
	refFsa <- system.file(package='SICtools','extdata','example.ref.fasta')
	
	## test deletion
	deletionTest <- indelDiff(bam1,bam2,refFsa,'chr04',962501,1026983,pValueCutOff=1,gtDistCutOff=0,verbose=FALSE)
	checkEquals(deletionTest$pos, 1026860)
	checkEquals(deletionTest$ref, 'TCG')
	checkEquals(deletionTest$altGt1, 'T')

	## test insertion
	insertionTest <- indelDiff(bam1,bam2,refFsa,'chr07',828514,828914,pValueCutOff=1,gtDistCutOff=0) 
	checkEquals(insertionTest$pos, 828714)
	checkEquals(insertionTest$ref, 'CAAAAAAA')
	checkEquals(insertionTest$altGt1, 'CAAAAAAAA')
	
	## test NULL
	checkTrue(is.null(indelDiff(bam1,bam2,refFsa,'chr04',962501,1026983,pValueCutOff=1,gtDistCutOff=2,verbose=FALSE)))
	
}




