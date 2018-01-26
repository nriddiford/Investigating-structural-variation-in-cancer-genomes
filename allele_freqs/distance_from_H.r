library(plotly)
library(ggplot2)
library(scales)

#setwd("/Users/Nick/iCloud/Desktop/snp_data/data")

setwd("/Users/Nick_curie/Desktop/script_test/snp_data")

path = getwd()

file.names <- dir(path, pattern =".txt")

for(i in 1:length(file.names)){
  cat("Processing", file.names[i])

  snps<-read.delim(file.names[i], header = T)

  parts<-strsplit(file.names[i], '[.]')[[1]]
  sample<-parts[1]
  
  p <- ggplot(data = snps, aes(x = position, y = freq, colour = type ), show.legend = FALSE)
  p <- p + stat_smooth(aes(fill = factor(type)), size = 0.25, alpha = 0.15, show.legend = FALSE)
  p <- p + geom_line()
  p <- p + scale_y_continuous(breaks = seq(0, 100, by = 25))
  p <- p + facet_wrap(~chrom, scale="free_x")
  p <- p + ggtitle( paste( sample ) )

  outfile <- paste(sample, '_', "allele_freqs", '.pdf', sep='')
  ggsave(outfile, scale = 0.9)
  
}

library(plotly)
library(ggplot2)
library(scales)

setwd("/Users/Nick_curie/Desktop/script_test/snp_data")

file<-'A512R21.dist_from_h.txt'

snps<-read.delim(file, header = T)
parts<-strsplit(file, '[.]')[[1]]
sample<-parts[1]
#
chrom2plot<-"2R"
#
chromosome<-subset(snps, snps$chrom  == chrom2plot )

# single chrom
# ggplot
p <- ggplot(data = chromosome, aes(x = position, y = dist_ratio), show.legend = FALSE)
p <- p + stat_smooth(aes(fill = dist_ratio), size = 0.25, alpha = 0.15, show.legend = FALSE)
p <- p + geom_point(aes(colour = -dist_ratio, text = paste("count: ", count )), show.legend = FALSE) 
# g <- g + scale_y_continuous(breaks = seq(0, 50, by = 10))
p <- p + ggtitle( paste( sample, chrom2plot ) )

p

# # plotly
g <- ggplot(data = chromosome, aes(x = position, y = dist_ratio), show.legend = FALSE)
g <- g + stat_smooth(aes(fill = dist_ratio), size = 0.25, alpha = 0.15, show.legend = FALSE)
g <- g + geom_point(aes(colour = -dist_ratio, text = paste("count: ", count )), show.legend = FALSE) 
# g <- g + scale_y_continuous(breaks = seq(0, 50, by = 10))
g <- g + ggtitle( paste( sample, chrom2plot ) )

ggplotly(g)

# all chroms

p <- ggplot(data = snps, aes(x = position, y = dist_ratio), show.legend = FALSE)
p <- p + stat_smooth(aes(fill = dist_ratio), size = 0.25, alpha = 0.15, show.legend = FALSE)
p <- p + geom_point(aes(colour = -dist_ratio, text = paste("count: ", count )), show.legend = FALSE) 
p <- p + facet_wrap(~chrom, scale="free_x")
p <- p + ggtitle( paste( sample ) )

p

 # plotly
g <- ggplot(data = snps, aes(x = position, y = freq, colour = type ), show.legend = FALSE)
g <- g + stat_smooth(aes(fill = factor(type)), size = 0.25, alpha = 0.15, show.legend = FALSE)
g <- g + geom_line()
g <- g + scale_y_continuous(breaks = seq(0, 100, by = 25))
g <- g + facet_wrap(~chrom, scale="free_x")
g <- g + ggtitle( paste( sample ) )

ggplotly(p)
p
p
