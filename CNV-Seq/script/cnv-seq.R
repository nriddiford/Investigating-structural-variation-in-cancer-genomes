
library(dplyr)

file.cleanR <- function(cnv_file) {

    clean_file <- filter(cnv_file, is.finite(log2) )
    
    cat("Filtering on number of reads per window (each window must have more than 10 reads in tum & cont)", "\n")
    clean_file <- filter(clean_file, test > 10 & ref > 10)
  
}

cnv.print <- function(cnv, file="")
{
	cat('cnv', 'chromosome', 'start', 'end', 'size', 'log2', 'p.value', sep="\t", file=file,fill=TRUE, append=FALSE)
	for(i in seq(max(min(cnv$cnv),1), max(cnv$cnv)))
	{
		sub <- subset(cnv, cnv==i)
		start <- ceiling(mean(c(min(sub$start), min(sub$position))))
		end <- floor(mean(c(max(sub$end), max(sub$position))))
		
		if(is.infinite(start)) next
                if(is.infinite(end)) next

		cat(paste('CNVR_',i,sep=''), paste('chr', unique(sub$chromosome), sep=''), start, end, end-start+1, unique(sub$cnv.log2), unique(sub$cnv.p.value), sep="\t", file=file, fill=TRUE, append=TRUE)
	}
}

cnv.summary <- function(cnv)
{
	true <- subset(cnv, cnv>0)
	
	max.size <- max(true$cnv.size)
	min.size <- min(true$cnv.size)
	count <- length(unique(true$cnv))
	nt.size <- 0
	for(i in unique(true$cnv)) nt.size <- nt.size + max(true[which(true$cnv==i),]$cnv.size)
	mean.size <- round(nt.size/count, 0)
	median.size <- median(unique(true$cnv.size))
	percen <- round(100*nrow(true)/nrow(cnv),1)
	
	cat('CNV percentage in genome: ', percen, "%\n", sep='')
	cat('CNV nucleotide content:', nt.size, fill=TRUE)
	cat('CNV count:', count, fill=TRUE)
	cat('Mean size:', mean.size, fill=TRUE)
	cat('Median size:', median.size, fill=TRUE)
	cat('Max Size:', max.size, fill=TRUE)
	cat('Min Size:', min.size, fill=TRUE)
}


args <- commandArgs()

cnv_file<-args[6]
id<-args[7]

cat("Reading .cnv file for ", id, ": ", cnv_file, "\n", sep='')


read_file_in<-read.delim(cnv_file)

clean_file <- file.cleanR(read_file_in)


cnv.print(clean_file, file=paste(id, '_cnvs.txt', sep=''))

