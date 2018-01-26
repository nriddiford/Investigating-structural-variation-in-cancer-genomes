library(plotly)
library(ggplot2)
library(scales)

setwd("/Users/Nick_curie/Desktop/script_test/snp_data/data")

# args <- commandArgs()
#
# file<-args[6]
#
# sample<- args[7]
# chrom2plot<- args[7]

file<-'A512R21.tagged.nodups.SC_q:15_cov:20_tp:0.5.snp'

snps<-read.delim(file, header = T)
parts<-strsplit(file, '[.]')[[1]]
sample<-parts[1]

chrom2plot<-"3R"

snps$normal_var_freq<- as.numeric(gsub("%", "", snps$normal_var_freq))
snps$tumor_var_freq<- as.numeric(gsub("%", "", snps$tumor_var_freq))

# roh <- with(snps, normal_var_freq > 95 & tumor_var_freq > 95)
# snps$somatic_status <- as.character(snps$somatic_status)
# snps$somatic_status[roh] <- "ROH"

chromosome<-subset(snps, snps$chrom  == chrom2plot )

allele_ratio<-((chromosome$tumor_var_freq)/(chromosome$normal_var_freq))

allele_ratio[is.infinite(allele_ratio)] <- 0

p <- ggplot(data = chromosome, aes(x = position, y = allele_ratio))

p <- p + geom_point(aes(colour = somatic_status), alpha=0.5)

p <- p + ggtitle( paste( sample, chrom2plot ) )

p <- p + geom_hline(yintercept=1, colour="slateblue", alpha=.7, linetype="dotted")
p


# plotly version
p <- ggplot(data = chromosome, aes(x = position, y = allele_ratio))

p <- p + geom_point(aes(colour = somatic_status,
                        alpha = 0.5,
                        text  = paste("somatic status: ", somatic_status,
                                      "<br>tum freq: "  , tumor_var_freq, "%",
                                      "<br>cont freq: " , normal_var_freq, "%")
                        )
                    )

p <- p + ggtitle( paste( sample, chrom2plot ) )
p <- p + geom_hline(yintercept=1, colour="slateblue", alpha=.7, linetype="dotted")

ggplotly(p)
