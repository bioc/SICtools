indelDiff <-
function(bam1,bam2,refFsa,regChr,regStart,regEnd,minBaseQuality = 13,
				minMapQuality = 0,nCores = 1,pValueCutOff = 0.05,gtDistCutOff = 0.1,verbose=TRUE){
	
	## set NULL
	rowIndex <- NULL
	
	## validate the parameters
	if(!file.exists(bam1)){stop("ERROR: input bam1 file doesn't exist!")}
	if(!file.exists(bam2)){stop("ERROR: input bam2 file doesn't exist!")}
	if(!file.exists(refFsa)){stop("ERROR: input refFsa file doesn't exist!")}
	
	if(regChr < 0){stop("ERROR: please input valid regChr!")}
	if(as.numeric(regStart) < 0){stop("ERROR: regStart should be more than 0!")}
	if(as.numeric(regEnd) < 0){stop("ERROR: regEnd should be more than 0!")}
	if(as.numeric(regStart) > as.numeric(regEnd)){stop("ERROR: regStart should be smaller than regEnd!")}
	
	regStart <- format(as.numeric(regStart),scientific=FALSE)
	regEnd <- format(as.numeric(regEnd),scientific=FALSE)
	minBaseQuality <- as.numeric(minBaseQuality)
	minMapQuality <- as.numeric(minMapQuality)
	pValueCutOff <- as.numeric(pValueCutOff)
	gtDistCutOff <- as.numeric(gtDistCutOff)
	
	## number of cores used
	registerDoMC(cores=nCores)
	
	## get mpileupPlus result; the samtools will calculate BAQ for SNP around indel
	pathSICtools <- system.file(package = "SICtools","etc","samtools2SIC")
	
	ppIndel <- read.delim(pipe(paste(pathSICtools,' mpileup -Q ',minBaseQuality,' -q ',minMapQuality, ' -Ogf ',refFsa,' ',bam1,' ',bam2,' -r ',regChr,':',regStart,'-',regEnd,sep='')),
			header=FALSE,colClasses = 'character',col.names=paste('V',1:12,sep=''))
	
	## remove items with too low indel hits
	ppIndel <- ppIndel[which(!(unlist(lapply(ppIndel$V6,
												function(x){sum(gregexpr('[+-]',x)[[1]] > 0)}))/as.numeric(ppIndel$V5) < 0.05 
								& unlist(lapply(ppIndel$V10,
												function(x){sum(gregexpr('[+-]',x)[[1]] > 0)}))/as.numeric(ppIndel$V9) < 0.05) 
							& as.numeric(ppIndel$V5) > 0 
							& as.numeric(ppIndel$V9) > 0),,drop=FALSE]
	
	if(nrow(ppIndel) > 0){ 
		
		## test out
		testOutIndel <- foreach(rowIndex=1:nrow(ppIndel)) %dopar% 
				.indelDiffFunc(ppIndel,rowIndex,refFsa,pValueCutOff,gtDistCutOff,verbose)
		
		## remove NULL
		testOutIndel <- testOutIndel[!sapply(testOutIndel,is.null)]
		
		## combine the result
		testOutIndel <- do.call(rbind,testOutIndel)
		
		if(!is.null(testOutIndel)){
			
			testIndelMtx <- testOutIndel[,-c(1,3:5),drop=FALSE]
			class(testIndelMtx) <- 'numeric'
			testIndelDf <- data.frame(chr=testOutIndel[,1],ref=testOutIndel[,3],alt1=testOutIndel[,4],alt2 = testOutIndel[,5],testIndelMtx,stringsAsFactors=FALSE)
			
			## at least one result
			if(nrow(testIndelDf) > 0){
				
				rownames(testIndelDf) <- 1:nrow(testIndelDf)
				testIndelDf <- testIndelDf[,c(1,5,2:4,6:13)]
				colnames(testIndelDf) <- c('chr','pos','ref','altGt1','altGt2','refBam1Count','altGt1Bam1Count',
						'altGt1Bam2Count','refBam2Count','altGt2Bam1Count','altGt2Bam2Count','p.value','d.value')
				
				testIndelDf <- testIndelDf[,c('chr','pos','ref','altGt1','altGt2','refBam1Count','altGt1Bam1Count','altGt2Bam1Count',
								'refBam2Count','altGt1Bam2Count','altGt2Bam2Count','p.value','d.value')]
				
				return(testIndelDf[,,drop=FALSE])
				
			}
		}
	}
}
