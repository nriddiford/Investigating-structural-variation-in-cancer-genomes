library(ggplot2)


snps<-read.table("HUM-7_somatic.data.txt", header = FALSE)
colnames(snps)=c("chroms","trans","freq")



### Barchart

# plot all on same row
#ggplot(snps, aes(x = trans, y = freq, group = chroms, fill = trans)) + geom_bar(position="dodge",stat="identity")+facet_grid(chroms~.)

# Displays grids by chromosome
#ggplot(snps, aes(x = trans, y = freq, group = chroms, fill = trans)) + geom_bar(position="dodge",stat="identity")+facet_wrap(~ chroms ) + theme(axis.text.x = element_text(angle=90))

# Displays grids by mutation

ggplot(snps, aes(x = chroms, y = freq, group = trans, fill = chroms)) + geom_bar(position="dodge",stat="identity")+facet_wrap(~ trans )

ggplot(snps, aes(x = chroms, y = freq, group = trans, fill = chroms)) + ylim(0, 50) + geom_bar(position="dodge",stat="identity")+facet_wrap(~ trans ) + theme(text = element_text(size=20))



### Boxplot

colnames(snps)=c("chroms","trans", "count", "freq")
ggplot(snps, aes(x = chroms, y = count, group = trans, fill = chroms)) + geom_dotplot(position="dodge",stat="identity")+facet_wrap(~ trans )


ggplot(snps, aes(x = chroms, y = count)) + geom_boxplot() +facet_wrap(~ trans )
