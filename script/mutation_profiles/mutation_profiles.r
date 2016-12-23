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






# For GW trinucs;
# Works for all data
ggplot(snps, aes(x = trans, y = freq, group = tri, fill = tri)) + geom_bar(position="dodge",stat="identity")+facet_wrap(~ trans)

# Try to get individual plots per transition
ggplot(snps, aes(x = tri, y = freq, group = trans, fill = trans)) + geom_bar(position="dodge",stat="identity") + theme(axis.text.x = element_text(angle=90))

ggplot(snps, aes(tri, freq)) + geom_bar(position="dodge",stat="identity")




c_t_2 <- read.table("C_T_muts.txt", header = FALSE)
colnames(c_t_2)=c("tri", "trans", "freq", "sample")
head(c_t_2)

 g <- ggplot(c_t_2, aes(tri,freq, group = sample)) + geom_point(aes(colour = sample))
g <- ggplot(c_t_2, aes(tri,freq, group = sample)) + geom_point(aes(colour = sample)) + geom_bar(position="dodge",stat="identity", alpha = 0.1)
g <- ggplot(c_t_2, aes(tri,freq, group = sample)) + geom_point(aes(colour = sample))+facet_wrap(~ trans) + theme(axis.text.x = element_text(angle=90))


Hum1 <- read.table("HUM-1.TEx.q4.frba.GWtrinucs.txt", header = FALSE)
colnames(Hum1)=c("tri", "trans", "freq", "sample")
head(Hum1)
ggplot(Hum1, aes(tri,freq, group = sample)) + geom_point(aes(colour = sample))+facet_wrap(~ trans,scale="free_x") + theme(axis.text.x = element_text(angle=90))

