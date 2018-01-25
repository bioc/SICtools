snpDiff <-
		function(bam1,bam2,refFsa,regChr,regStart,regEnd,minBaseQuality = 13,
				minMapQuality = 0,nCores = 1,pValueCutOff= 0.05,baseDistCutOff = 0.1,verbose=TRUE){
	
	## disable dimention dropping for matrix
	`[` <- function(...) base::`[`(...,drop=FALSE)
	
	## validate the parameters
	if(!file.exists(bam1)){stop("ERROR: input bam1 file doesn't exist!")}
	if(!file.exists(bam2)){stop("ERROR: input bam2 file doesn't exist!")}
	if(!file.exists(refFsa)){stop("ERROR: input refFsa file doesn't exist!")}
	
	if(regChr < 0 ){stop("ERROR: please input valid regChr!")}
	regStart <- as.numeric(regStart);if(regStart < 0 ){stop("ERROR: regStart should be more than 0!")}
	regEnd <- as.numeric(regEnd);if(regEnd < 0 ){stop("ERROR: regEnd should be more than 0!")}
	if(regStart > regEnd){stop("ERROR: regStart should be smaller than regEnd!")}
	
	pValueCutOff <- as.numeric(pValueCutOff)
	baseDistCutOff <- as.numeric(baseDistCutOff) 
	minBaseQuality <- as.numeric(minBaseQuality)
	minMapQuality <- as.numeric(minMapQuality)
	
	## number of cores used
	registerDoParallel(cores=nCores)
	
	################################ function to call difference ##############################
	calcInfoRange <- function(x){
		
		countMtx <- matrix(x[['seq']],ncol=10,byrow=TRUE)
		countMtx <- cbind(countMtx,x[['pos']])
		
		## duplicated index
		dupIndex <- duplicated(countMtx[,-c(5,10,11)])
		
		## max bam1 == max bam2 index
		maxBaseIndex <- rowMaxs(countMtx[,1:4])/rowSums(countMtx[,1:4]) > 0.95 & 
				rowMaxs(countMtx[,6:9])/rowSums(countMtx[,6:9]) > 0.95 & 
				max.col(countMtx[,1:4]) == max.col(countMtx[,6:9])
		
		## all 0 test
		allZeroIndex <- rowSums(countMtx[,1:4]) == 0 | rowSums(countMtx[,6:9]) == 0		
		
		## test index
		testIndex <- !dupIndex & !maxBaseIndex & !allZeroIndex
		
		if(sum(testIndex) > 0){ ## at least one candidate
			
			countMtxTest <- countMtx[testIndex,]
			
			testResultMtx <- do.call(rbind,llply(1:nrow(countMtxTest),function(rowIndex){
						
						rowInfo <- countMtxTest[rowIndex,]
						
						## get test matrix
						testMtx <- matrix(rowInfo[c(1:4,6:9)],nrow=2,byrow=TRUE)
						
						## remove 0 columns
						testMtx <- testMtx[,colSums(testMtx) > 0]
						
						## get pValue and dist
						pValue <- fisher.test(testMtx)$p.value
						baseDist <- as.numeric(dist(testMtx/rowSums(testMtx),method = 'euclidean'))
						
						pValue <- ifelse(is.na(pValue),1,pValue)
						baseDist <- ifelse(is.na(baseDist),0,baseDist)
						
						if(pValue <= as.numeric(pValueCutOff) & baseDist >= as.numeric(baseDistCutOff)){
							
							c(rowInfo[1],rowInfo[2],rowInfo[3],rowInfo[4],rowInfo[5],
											rowInfo[6],rowInfo[7],rowInfo[8],rowInfo[9],rowInfo[10],rowInfo[11],
											pValue,baseDist)
						}
					}))
						
			## for genotype duplicated positions, trace back
			countMtxDup <- countMtx[!maxBaseIndex & !allZeroIndex,]
			countMtxDupIndex <- duplicated(countMtxDup[,-11])
			
			if(sum(countMtxDupIndex) > 0 & !is.null(testResultMtx)){
				
				countMtxDupFull <- cbind(countMtxDup,testResultMtx[match(apply(countMtxDup[,1:10],1,paste,collapse='_'),apply(testResultMtx[,1:10],1,paste,collapse='_')),c(12,13)])
				countMtxDupFull <- countMtxDupFull[!is.na(countMtxDupFull[,12]),]
				testResultMtx <- unique(rbind(testResultMtx,countMtxDupFull))
			}
			
			return(testResultMtx)          
		}
	}
	
	############################################################################
	
	
	## build bam files for pileup
	ppFiles <- PileupFiles(c(bam1,bam2))
	
	## split into 10M block
	regSplit <- shift(as(breakInChunks(regEnd - regStart + 1, chunksize=1e7L),'IRanges'),regStart-1)
	
	## for each region
	regSplitDiffFunc <- function(splitIndex){
		
		## print progress
		if(verbose){
			print(paste(regChr,' ',start(regSplit)[splitIndex],' ',end(regSplit)[splitIndex],sep=''))
		}
		
		regSplitParam <- ApplyPileupsParam(flag = scanBamFlag(isDuplicate=FALSE),
				which=GRanges(regChr, IRanges(start=start(regSplit[splitIndex]), end=end(regSplit[splitIndex]))), 
				what='seq',minDepth = 0L,minBaseQuality, minMapQuality,maxDepth = 2500L, yieldBy = 'range')
		
		applyPileups(ppFiles,calcInfoRange,param=regSplitParam)[[1]]
		
	}
	
	## for R CMD check
	splitIndex <- NULL
	
	## range out
	rangeOut <- ldply(1:length(regSplit),regSplitDiffFunc,.parallel=TRUE)
		
		if(nrow(rangeOut) > 0){
			
			colnames(rangeOut) <- c('A1','C1','G1','T1','N1','A2','C2','G2','T2','N2','pos','p.value','d.value')
			rangeOut$chr <- regChr
			
			## get ref by Biostrings
			rangeOut$ref <- as.character(getSeq(FaFile(refFsa),GRanges(seqnames=rangeOut$chr,IRanges(start=rangeOut$pos,width=1))))
			rangeOut <- rangeOut[order(rangeOut$pos),]
			
			## format output    
			rangeOut <- rangeOut[,c('chr','pos','ref','A1','C1','G1','T1','N1','A2','C2','G2','T2','N2','p.value','d.value')]
			return(rangeOut)
		}else{
			NULL
		}
}
