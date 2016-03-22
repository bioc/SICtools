snpDiff <-
function(bam1,bam2,refFsa,regChr,regStart,regEnd,minBaseQuality = 13,
		minMapQuality = 0,nCores = 1,pValueCutOff= 0.05,baseDistCutOff = 0.1,verbose=TRUE){
	
	## validate the parameters
	if(!file.exists(bam1)){stop("ERROR: input bam1 file doesn't exist!")}
	if(!file.exists(bam2)){stop("ERROR: input bam2 file doesn't exist!")}
	if(!file.exists(refFsa)){stop("ERROR: input refFsa file doesn't exist!")}
	
	if(regChr < 0 ){stop("ERROR: please input valid regChr!")}
	if(as.numeric(regStart) < 0 ){stop("ERROR: regStart should be more than 0!")}
	if(as.numeric(regEnd) < 0 ){stop("ERROR: regEnd should be more than 0!")}
	if(as.numeric(regStart) > as.numeric(regEnd)){stop("ERROR: regStart should be smaller than regEnd!")}
	
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
		dupIndex <- duplicated(countMtx[,-c(5,10,11),drop=FALSE])
		
		## max bam1 == max bam2 index
		maxBaseIndex <- rowMaxs(countMtx[,1:4,drop=FALSE])/rowSums(countMtx[,1:4,drop=FALSE]) > 0.95 & 
				rowMaxs(countMtx[,6:9,drop=FALSE])/rowSums(countMtx[,6:9,drop=FALSE]) > 0.95 & 
				max.col(countMtx[,1:4,drop=FALSE]) == max.col(countMtx[,6:9,drop=FALSE])
		
#		maxBaseIndex <- rowMaxs(countMtx[,,drop=FALSE],cols=1:4)/rowSums(countMtx[,,drop=FALSE],cols=1:4) > 0.95 & 
#				rowMaxs(countMtx[,,drop=FALSE],cols=6:9)/rowSums(countMtx[,,drop=FALSE],cols=6:9) > 0.95 & 
#				colMaxs(countMtx[,,drop=FALSE],cols=1:4) == colMaxs(countMtx[,,drop=FALSE],cols=6:9)
		
		## all 0 test
		allZeroIndex <- rowSums(countMtx[,1:4,drop=FALSE]) == 0 | rowSums(countMtx[,6:9,drop=FALSE]) == 0
#		allZeroIndex <- rowSums(countMtx[,,drop=FALSE],cols=1:4) == 0 | rowSums(countMtx[,,drop=FALSE],cols=6:9) == 0
		
		
		## test index
		testIndex <- !dupIndex & !maxBaseIndex & !allZeroIndex
		
		if(sum(testIndex) > 0){ ## at least one candidate
			
			countMtxTest <- countMtx[testIndex,,drop=FALSE]
			
			testTmp <- lapply(1:nrow(countMtxTest),function(rowIndex){
						
						rowInfo <- countMtxTest[rowIndex,]
						
						## get test matrix
						testMtx <- matrix(rowInfo[c(1:4,6:9)],nrow=2,byrow=TRUE)
						
						## remove 0 columns
						testMtx <- testMtx[,colSums(testMtx) > 0,drop=FALSE]
						
						## get pValue and dist
						pValue <- fisher.test(testMtx)$p.value
						baseDist <- as.numeric(dist(testMtx/rowSums(testMtx),method = 'euclidean'))
						
						pValue <- ifelse(is.na(pValue),1,pValue)
						baseDist <- ifelse(is.na(baseDist),0,baseDist)
						
						if(pValue <= as.numeric(pValueCutOff) & baseDist >= as.numeric(baseDistCutOff)){
							
							as.character(c(names(x[['seqnames']]),rowInfo[11],
											rowInfo[1],rowInfo[2],rowInfo[3],rowInfo[4],rowInfo[5],
											rowInfo[6],rowInfo[7],rowInfo[8],rowInfo[9],rowInfo[10],
											format(pValue,digits=3,scientific=TRUE),
											format(baseDist,digits=3,scientific=TRUE)))
						}
					})
			
			testTmp2 <- do.call(rbind,testTmp[,drop=FALSE])
			testMtx <- testTmp2[,-1,drop=FALSE]
			class(testMtx) <- 'numeric'
			testDf <- data.frame(chr = testTmp2[,1],testMtx,stringsAsFactors=FALSE)
			
			## for genotype duplicated positions, trace back
			countMtxDup <- countMtx[!maxBaseIndex & !allZeroIndex,,drop=FALSE]
			countMtxDupDf <- as.data.frame(countMtxDup[duplicated(countMtxDup[,-11,drop=FALSE]),,drop=FALSE],stringsAsFactors=FALSE)
			
			if(nrow(countMtxDupDf) > 0 & nrow(testDf) > 0){
				
				countMtxDupDfFull <- cbind(countMtxDupDf,testDf[match(apply(countMtxDupDf[,-11,drop=FALSE],1,paste,collapse=' '),apply(testDf[,3:12],1,paste,collapse=' ')),c('X12','X13')])
				countMtxDupDfFull <- countMtxDupDfFull[!is.na(countMtxDupDfFull$X12),,drop=FALSE]
				
				if(nrow(countMtxDupDfFull) > 0){
					countMtxDupDfFull[,'X14'] <- unique(testDf$chr)
					countMtxDupDfFull <- countMtxDupDfFull[,c(14,11,1:10,12,13)]
					colnames(countMtxDupDfFull) <- colnames(testDf)
					testDf <- rbind(testDf,countMtxDupDfFull)
					testDf <- testDf[order(testDf[,'X1']),]
				}
				
			}
			
			return(testDf)          
		}
	}
	
	############################################################################
	
	
	## build bam files for pileup
	ppFiles <- PileupFiles(c(bam1,bam2))
	
	## get region info
	regStart <- as.numeric(regStart)
	regEnd <- as.numeric(regEnd) + 1
	
	## split into 10M block
	regSplit <- shift(as(breakInChunks(regEnd - regStart,1e7L),'IRanges'),regStart-1)
	
	## for each region
	regSplitDiffFunc <- function(splitIndex){
		
		## print progress
		if(verbose){
			print(paste(regChr,' ',start(regSplit)[splitIndex],' ',end(regSplit)[splitIndex],sep=''))
		}
		
		regSplitParam <- ApplyPileupsParam(flag = scanBamFlag(isDuplicate=FALSE),
				which=GRanges(regChr, IRanges(start=start(regSplit)[splitIndex], end=end(regSplit)[splitIndex])), 
				what='seq',minDepth = 0L,minBaseQuality, minMapQuality,maxDepth = 2500L, yieldBy = 'range')
		
		applyPileups(ppFiles,calcInfoRange,param=regSplitParam)[[1]]
		
	}
	
	## for R CMD check
	splitIndex <- NULL
	
	## range out
	rangeOut <- foreach(splitIndex=1:length(regSplit)) %dopar% regSplitDiffFunc(splitIndex)
	
	## remove NULL
	rangeOut <- rangeOut[!sapply(rangeOut,is.null)]
	
	## combine the results
	rangeOut <- do.call(rbind,rangeOut)
	
	if(!is.null(rangeOut)){
		if(nrow(rangeOut) > 0){
			
			## get ref by Biostrings
			rangeOut$ref <- as.character(getSeq(FaFile(refFsa),GRanges(seqnames=rangeOut[,1],IRanges(start=rangeOut[,2],width=1))))
			
			## format output    
			rangeOut <- rangeOut[,c(1,2,15,3:14),drop=FALSE]
			colnames(rangeOut) <- c('chr','pos','ref','A1','C1','G1','T1','N1','A2','C2','G2','T2','N2','p.value','d.value')
			return(rangeOut[,,drop=FALSE])
		}}
}
