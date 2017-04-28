library(plotly)
library(ggplot2)
library(scales)
library(RColorBrewer)

file.clean <- function(cnv_file)
{
   	# Removes values of '-Inf'
 	clean_file<-cnv_file[!(cnv_file$log2=="-Inf"),]
	clean_file<-cnv_file[!(cnv_file$log2=="Inf"),]

	# Remove if fewer than 20 reads in window
  	clean_file<-cnv_file[!(cnv_file$test< 20),]
	clean_file<-cnv_file[!(cnv_file$ref< 20),]
}


######################
## All chroms batch ##
######################
 
plot.all.grid <- function(path=NA)
{
    if(is.na(path)) {
		path <- "cnvs/"
    	cat("No file path provided - defaulting to", path, "\n")
	}
	
	file.names <- dir(path, pattern =".cnv")
	
	for(i in 1:length(file.names)){
    	
		cat("Processing file", file.names[i], "\n")
		
  	  	parts<-strsplit(file.names[i], '[.]')[[1]]
  	  	sample<-parts[1]

  	  	read_file_in<-read.delim(paste("cnvs/", file.names[i], sep = ''), header = T)

		clean_file<-file.clean(read_file_in)
		
  	  	cols <- brewer.pal(n = 5, name = "RdBu") 
  
  	  	p <- ggplot()
  	  	p <- p + geom_point(data=clean_file, aes(x = start, y = log2, colour = log2), size = 1)
  	  	p <- p + ylim(-5,5)
  	  	p <- p + scale_colour_gradientn(colours = cols, 
	  	  	values = rescale(c(-2, -0.25, 0, 0.25, 2)),
      		guide = "colorbar", limits=c(-5, 5))
  		p <- p + facet_wrap(~chromosome, scale="free_x")
  		p <- p + ggtitle( paste( sample ) )
  	  	p <- p + theme(plot.title = element_text(hjust = 0.5))
  
  	  	outfile <- paste(sample, '_', "CNVs", '.pdf', sep='')
  
  	  	cat("Writing file", outfile, "to `../plots/`", "\n")
  	  	ggsave(paste("plots/", outfile, sep = ''), width = 20, height = 10)
  
	}
}

##################
## Single chrom ##
##################

plot.chrom <- function(chrom=NA, cnv_file)
{
  cat("Processing", cnv_file, "\n")
  
  if(is.na(chrom)) {
	  chrom <- "3R"
  	  cat("No chromosome specified - defaulting to", chrom, "\n")
  }
  else {
  	  cat("Plotting", chrom, "chromosome", "\n")
  }	  

	base=basename(cnv_file)
	parts<-strsplit(base, '[.]')[[1]]
	sample<-parts[1]

	read_file_in<-read.delim(cnv_file, header = T)

	clean_file<-file.clean(read_file_in)

	cols <- brewer.pal(n = 5, name = "RdBu") 
		
	chrom_data<-subset(clean_file, clean_file$chromosome == chrom)

	p <- ggplot()
	p <- p + geom_point(data=chrom_data, aes(x = start, y = log2, colour = log2))
	p <- p + scale_alpha(range = c(0.1, 5))
	p <- p + ylim(-5,5)
	
	p <- p + scale_colour_gradientn(colours = cols, 
	                                values = rescale(c(-2, -0.2, 0, 0.2, 2)),
	                                guide = "colorbar", limits=c(-5, 5))
	p <- p + ggtitle( paste( sample, " ", chrom, sep = '' ) )
	p <- p + theme(plot.title = element_text(hjust = 0.5))
	
	p <- p + geom_hline(yintercept=-0.322, colour="slateblue", alpha=.25, linetype = "dotted")
	p <- p + geom_hline(yintercept=0.322, colour="slateblue", alpha=.25, linetype = "dotted")
	p <- p + geom_hline(yintercept=-0.585, colour="blue", alpha=.3)
	p <- p + geom_hline(yintercept=0.585, colour="blue", alpha=.3)
	p <- p + geom_hline(yintercept=-1, colour="royalblue", alpha=.4)
	p <- p + geom_hline(yintercept=1, colour="royalblue", alpha=.4)
	
	outfile <- paste(sample, '_', chrom, '_', "CNVs", '.pdf', sep='')
	
	cat("Writing file", outfile, "to `../plots/`", "\n")
  	ggsave(paste("plots/", outfile, sep = ''), width = 20, height = 10)
	
}