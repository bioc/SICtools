indelDiff <-
function(bam1,bam2,refFsa,regChr,regStart,regEnd,minBaseQuality = 13,
				minMapQuality = 0,nCores = 1,pValueCutOff = 0.05,gtDistCutOff = 0.1,verbose=TRUE){
	
		
	## disable dimention dropping for matrix
	`[` <- function(...) base::`[`(...,drop=FALSE)
	
	## set NULL
	rowIndex <- NULL
	
	## validate the parameters
	if(!file.exists(bam1)){stop("ERROR: input bam1 file doesn't exist!")}
	if(!file.exists(bam2)){stop("ERROR: input bam2 file doesn't exist!")}
	if(!file.exists(refFsa)){stop("ERROR: input refFsa file doesn't exist!")}
	
	if(regChr < 0){stop("ERROR: please input valid regChr!")}
	regStart <- as.numeric(regStart);if(regStart < 0 ){stop("ERROR: regStart should be more than 0!")}
	regEnd <- as.numeric(regEnd);if(regEnd < 0 ){stop("ERROR: regEnd should be more than 0!")}
	if(regStart > regEnd){stop("ERROR: regStart should be smaller than regEnd!")}
	
	regStart <- format(regStart,scientific=FALSE)
	regEnd <- format(regEnd,scientific=FALSE)
	minBaseQuality <- as.numeric(minBaseQuality)
	minMapQuality <- as.numeric(minMapQuality)
	pValueCutOff <- as.numeric(pValueCutOff)
	gtDistCutOff <- as.numeric(gtDistCutOff)
	
	## number of cores used
	registerDoParallel(cores=nCores)
	
	## get mpileupPlus result; the samtools will calculate BAQ for SNP around indel
	pathSICtools <- system.file(package = "SICtools","etc","samtools2SIC")
	
	ppIndel <- read.delim(pipe(paste(pathSICtools,' mpileup -Q ',minBaseQuality,' -q ',minMapQuality, ' -Ogf ',refFsa,' ',bam1,' ',bam2,' -r ',regChr,':',regStart,'-',regEnd,sep='')),
			header=FALSE,colClasses = 'character',col.names=paste('V',1:12,sep=''))
	
	## remove items with too low indel hits	
	ppIndel <- ppIndel[which(!(unlist(lapply(gregexpr('[+-]',ppIndel$V6),function(x){sum(x > 0)}),use.names=FALSE)/as.numeric(ppIndel$V5) < 0.05 & 
								
	                           unlist(lapply(gregexpr('[+-]',ppIndel$V10),function(x){sum(x > 0)}),use.names=FALSE)/as.numeric(ppIndel$V9) < 0.05) & 
						       as.numeric(ppIndel$V5) > 0 & as.numeric(ppIndel$V9) > 0),]
	
	if(nrow(ppIndel) > 0){
		
		## test out
		testOutIndel <- ldply(1:nrow(ppIndel),function(rowIndex){.indelDiffFunc(ppIndel,rowIndex,refFsa,pValueCutOff,gtDistCutOff,verbose)},.parallel=TRUE)
		
		if(nrow(testOutIndel) > 0){
			
			colnames(testOutIndel) <- c('chr','pos','ref','altGt1','altGt2','refBam1Count','altGt1Bam1Count',
					'altGt2Bam1Count','refBam2Count','altGt1Bam2Count','altGt2Bam2Count','p.value','d.value')
			class(testOutIndel$pos) <- 'numeric'
			class(testOutIndel$refBam1Count) <- 'numeric'
			class(testOutIndel$altGt1Bam1Count) <- 'numeric'
			class(testOutIndel$altGt2Bam1Count) <- 'numeric'
			class(testOutIndel$refBam2Count) <- 'numeric'
			class(testOutIndel$altGt1Bam2Count) <- 'numeric'
			class(testOutIndel$altGt2Bam2Count) <- 'numeric'
			class(testOutIndel$p.value) <- 'numeric'
			class(testOutIndel$d.value) <- 'numeric'
			
			return(testOutIndel)
		}
		
	}else{
		NULL
	}
}
