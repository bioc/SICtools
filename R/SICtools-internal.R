.indelDiffFunc <-
function(ppIndel,rowIndex,refFsa,pValueCutOff,gtDistCutOff,verbose){
	
	## expression to evaluate
	str_eval <- function(x){
		eval(parse(text=x))
	}
	
	## get row data
	rowDat <- as.matrix(ppIndel[rowIndex,])
	
	## tips
	if(verbose){
		print(as.character(rowDat[,2:4]))
	}
	
	## indel region
	indelLen <- max(as.numeric(strsplit(rowDat[1],',')[[1]]))
	
	## original genotypes of two bam files
	gtBam1 <- toupper(strsplit(rowDat[6],'\\|')[[1]])
	gtBam1[gtBam1 == '.'|gtBam1 == ','] <- '*'
	
	gtBam2 <- toupper(strsplit(rowDat[10],'\\|')[[1]])
	gtBam2[gtBam2 == '.'|gtBam2 == ','] <- '*'
	
	## indel length 
	indelLenBam1 <- str_extract(gtBam1,'[+-]\\d+')
	indelLenBam1[is.na(indelLenBam1)] <- '+0'
	
	indelLenBam2 <- str_extract(gtBam2,'[+-]\\d+')
	indelLenBam2[is.na(indelLenBam2)] <- '+0'
	
	## first parse soft clip on both sides
	gtBam1SoftLeft <- as.numeric(sub('S','',str_extract(strsplit(rowDat[8],',')[[1]][-1],'^\\d+S')))
	gtBam1SoftLeft[is.na(gtBam1SoftLeft)] <- 0
	
	gtBam1qLen <- sapply(str_extract_all(strsplit(rowDat[8],',')[[1]][-1],'\\d+(I|M)'),
			function(x){sum(as.numeric(gsub('I|M','',x)))})
	
	gtBam2SoftLeft <- as.numeric(sub('S','',str_extract(strsplit(rowDat[12],',')[[1]][-1],'^\\d+S')))
	gtBam2SoftLeft[is.na(gtBam2SoftLeft)] <- 0
	
	gtBam2qLen <- sapply(str_extract_all(strsplit(rowDat[12],',')[[1]][-1],'\\d+(I|M)'),
			function(x){sum(as.numeric(gsub('I|M','',x)))})
	
	## reads should cover more than 2bp of indel region based on insertion or deletion
	qposBam1 <- as.numeric(strsplit(rowDat[7],',')[[1]][-1])
	bam1Sel <- (qposBam1 - gtBam1SoftLeft >= 2) & 
			(unlist(lapply(paste((qposBam1 - gtBam1SoftLeft) + indelLen + 2, indelLenBam1,sep=''),str_eval)) <= gtBam1qLen)
	gtBam1 <- gtBam1[bam1Sel]
	gtBam1 <- sub('^[.,]','',gtBam1)
	
	qposBam2 <- as.numeric(strsplit(rowDat[11],',')[[1]][-1])
	bam2Sel <- (qposBam2 - gtBam2SoftLeft >= 2) & 
			(unlist(lapply(paste((qposBam2 - gtBam2SoftLeft) + indelLen + 2, indelLenBam2,sep=''),str_eval)) <= gtBam2qLen)
	gtBam2 <- gtBam2[bam2Sel]
	gtBam2 <- sub('^[.,]','',gtBam2)
	
	## remove same gtBam and at least one result
	gtSameIndex <- FALSE
	if(length(unique(gtBam1)) == 1 & length(unique(gtBam2)) == 1){gtSameIndex <- unique(gtBam1) == unique(gtBam2)}
	
	if(length(gtBam1) > 0 & length(gtBam2) > 0 & !gtSameIndex){
		
		## genotyping
		gtMtx <- rbind(table(factor(gtBam1,levels=unique(c(gtBam1,gtBam2)))),
				table(factor(gtBam2,levels=unique(c(gtBam1,gtBam2)))))
		
		## fill the reference allele if not exist
		if(sum(colnames(gtMtx) == '*') == 0) {gtMtx <- cbind(gtMtx,'*'=c(0,0))}
		
		## select two most alts
		if(ncol(gtMtx) >= 3){  
			
			gtCountMtx <- gtMtx[,c('*',names(sort(colSums(gtMtx[,which(colnames(gtMtx)  != '*')]),
											decreasing=TRUE)[1:2]))]     
			
		}else if(ncol(gtMtx) == 2){
			
			gtCountMtx <- gtMtx
			gtCountMtx <- gtCountMtx[,sort(colnames(gtCountMtx))]
			
		}
		
		## test for p.value and distance
		gtPvalue <- fisher.test(gtCountMtx)$p.value
		gtDist <- as.numeric(dist(gtCountMtx/rowSums(gtCountMtx),method = 'euclidean'))
		gtPvalue <- ifelse(is.na(gtPvalue),1,gtPvalue)
		gtDist <- ifelse(is.na(gtDist),0,gtDist)
		
		if(gtPvalue <= pValueCutOff & gtDist >= gtDistCutOff){
			
			## format genotype to output
			delSeqLen <- as.numeric(str_extract(colnames(gtCountMtx),'-\\d+'))
			delSeqLen[is.na(delSeqLen)] <- 0
			colnames(gtCountMtx)[delSeqLen < 0] <- '*'
			
			gtOut <- paste(rowDat[4],colnames(gtCountMtx),
					getSeq(FaFile(refFsa),GRanges(seqnames=rowDat[2],IRanges(as.numeric(rowDat[3])+1,
											width=indelLen+delSeqLen))), sep='')
			
			indelSNPIndex <- grepl('^[ATCGN]',colnames(gtCountMtx))
			
			## sometimes there will be SNP around indel
			if(sum(indelSNPIndex) > 0){
				
				gtOut[indelSNPIndex] <- colnames(gtCountMtx)[indelSNPIndex]
				
			}
			colnames(gtCountMtx) <- gsub('\\*|\\+\\d+','',gtOut)
			
			if(ncol(gtCountMtx) == 3){
				
				## output table
				tmpOut <- c(rowDat[,2:3],colnames(gtCountMtx),
								gtCountMtx[1,1], gtCountMtx[1,2],gtCountMtx[1,3],
								gtCountMtx[2,1], gtCountMtx[2,2],gtCountMtx[2,3],
								pValue = gtPvalue,
								gtDist = gtDist)
						
			}else if(ncol(gtCountMtx) < 3){
				
				## output table
				tmpOut <- c(rowDat[,2:3],colnames(gtCountMtx),NA,
								gtCountMtx[1,1], gtCountMtx[1,2],NA,
								gtCountMtx[2,1], gtCountMtx[2,2],NA,
								pValue = gtPvalue,
								gtDist = gtDist)
			}
			
			return(tmpOut)
		}
	}
}
