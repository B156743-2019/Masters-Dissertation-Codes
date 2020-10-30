#make accumulation plot
accm<- function () {
sp1 <- specaccum(t(otu_table(physeq1)))
#sp2 <- specaccum(t(otu_table(physeq1)), "random")
pdf('accumulation_plot2.pdf')
plot(sp1, ci=2, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
#boxplot(sp2, col="yellow", add=TRUE, pch="+")
dev.off()
}
