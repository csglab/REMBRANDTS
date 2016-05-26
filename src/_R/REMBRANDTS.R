library(gplots)

args <- commandArgs(trailingOnly = TRUE)
jobID <- args[1]
stringency <- as.numeric(args[2])

tempFolder <- paste("./tmp/",jobID,sep="")
outputFolder <- paste("./out/",jobID,sep="")

# read the exon and intron centered normalized reads, and the total read counts
exon <- read.csv(paste(tempFolder,"/vsd_normalized.exonic.all.centered.mx.txt",sep=""),sep="\t")
nSample <- ncol(exon)-1
exonRawCounts <- read.csv(paste(tempFolder,"/vsd_normalized.exonic.all.mx.txt",sep=""),sep="\t")
if( nSample != ncol(exonRawCounts)-1 )
  print("Error: unequal sample numbers")
intron <- read.csv(paste(tempFolder,"/vsd_normalized.intronic.all.centered.mx.txt",sep=""),sep="\t")
if( nSample != ncol(intron)-1 )
  print("Error: unequal sample numbers")
intronRawCounts <- read.csv(paste(tempFolder,"/vsd_normalized.exonic.all.mx.txt",sep=""),sep="\t")
if( nSample != ncol(intronRawCounts)-1 )
  print("Error: unequal sample numbers")

# calculate the median counts
exonRawCounts$MedianExon <- apply(exonRawCounts[,2:ncol(exonRawCounts)],1,median)
intronRawCounts$MedianIntron <- apply(intronRawCounts[,2:ncol(intronRawCounts)],1,median)

merged <- merge(exon,intron,by="GeneID")
merged <- merge( merged, exonRawCounts[,c(1,ncol(exonRawCounts))],by="GeneID")
merged <- merge( merged, intronRawCounts[,c(1,ncol(intronRawCounts))],by="GeneID")

# optimize the total read count threshold so as to maximize correlation between exon and intron fold-changes
for( i in seq(0,50,by=0.1))
	if( sum(merged$MedianIntron > i & merged$MedianExon > i) > 500 )
	{	
		correl <- cor(
			unlist(merged[ merged$MedianIntron > i & merged$MedianExon > i , 2:(nSample+1)]),
			unlist(merged[ merged$MedianIntron > i & merged$MedianExon > i, (nSample+2):(2*nSample+1)]) )

		if( i == 0 )
			correl_max <- correl_all <- correl
		
		if( correl_max < correl )
			correl_max <- correl
		
		if( correl >= (correl_max-correl_all)*stringency + correl_all )
			threshold <- i

		print( paste( i, correl, sum( merged$MedianIntron > i & merged$MedianExon > i) ) )
  
	}

print( paste( "Selected threshold is ", threshold, sep="" ) )

# the filtered exon counts
y <- unlist(merged[ merged$MedianIntron > threshold & merged$MedianExon > threshold, 2:(nSample+1)])
# the filtered intron counts
x <- unlist(merged[ merged$MedianIntron > threshold & merged$MedianExon > threshold, (nSample+2):(2*nSample+1)])

jpeg(file=paste(outputFolder,"/scatterplot.jpg",sep=""),
	width=300,height=350)
smoothScatter( x, y, colramp=colorRampPalette(c( "white", "red", "black" ) ), nbin=400, nrpoints=0, xlab="Log2 fold-change of intronic reads", ylab="Log2 fold-change of exonic reads" )

# The intron counts represent the effect of transcriptional changes
# The exon counts represent the combined effect of transcriptional and post-transcriptional changes
# Model the exon counts as a function of intron counts, so that by subtracting the fitted values
#  from exon counts, only the effect of post-transcriptional regulation remains

# create a copy of exon counts
exon.counts <- merged[ merged$MedianIntron > threshold & merged$MedianExon > threshold, 1:(nSample+1)]
# plot the similarity heatmap for exon counts
sim <- cor(as.matrix(exon.counts[,2:(nSample+1)]))
width  <- ncol(sim) * 36.9 + 231
height <- nrow(sim) * 36.9 + 231
jpeg(file=paste(outputFolder,"/exonic.filtered.correl.heatmap.jpg",sep=""),
	width=width,height=height)
heatmap.2(1-sim, distfun=as.dist, trace="none", margin=c(13, 13),symm=T,revC=T,breaks=seq(0,2,length.out=256),col=colorRampPalette(c("red","white","blue"))(255))

# create a copy of intron counts
intron.counts <- merged[ merged$MedianIntron > threshold & merged$MedianExon > threshold, c(1,(nSample+2):(2*nSample+1))]
# plot the similarity heatmap for intron counts
sim <- cor(as.matrix(intron.counts[,2:(nSample+1)]))
width  <- ncol(sim) * 36.9 + 231
height <- nrow(sim) * 36.9 + 231
jpeg(file=paste(outputFolder,"/intronic.filtered.correl.heatmap.jpg",sep=""),
	width=width,height=height)
heatmap.2(1-sim, distfun=as.dist, trace="none", margin=c(13, 13),symm=T,revC=T,breaks=seq(0,2,length.out=256),col=colorRampPalette(c("red","white","blue"))(255))

# sample-by-sample, deconvolute the PTR effect
ptr <- exon.counts
for( i in 2:(nSample+1) )
{
  fit <- glm( exon.counts[,i] ~ intron.counts[,i] )
  ptr[,i] = exon.counts[,i] - ( fit$coefficients[1] + fit$coefficients[2] * intron.counts[,i] )
}

# plot the similarity heatmap for deconvoluted PTR matrix
sim <- cor(as.matrix(ptr[,2:(nSample+1)]))
width  <- ncol(sim) * 36.9 + 231
height <- nrow(sim) * 36.9 + 231
jpeg(file=paste(outputFolder,"/stability.filtered.correl.heatmap.jpg",sep=""),
	width=width,height=height)
heatmap.2(1-sim, distfun=as.dist, trace="none", margin=c(13, 13),symm=T,revC=T,breaks=seq(0,2,length.out=256),col=colorRampPalette(c("red","white","blue"))(255))

# write the tables
write.table(exon.counts,paste(outputFolder,"/exonic.filtered.mx.txt",sep=""),sep="\t",quote=F,row.names=F)
write.table(intron.counts,paste(outputFolder,"/intronic.filtered.mx.txt",sep=""),sep="\t",quote=F,row.names=F)
write.table(ptr,paste(outputFolder,"/stability.filtered.mx.txt",sep=""),sep="\t",quote=F,row.names=F)
