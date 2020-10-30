#remove the fake sample '1000'
physeq<-phyloseq(
otu_table(physeq1)[,colnames(otu_table(physeq1))!='1000'],
sample_data(physeq1)[rownames(sample_data(physeq1))!='1000'],
tax_table(physeq1))
#filter samples
qp <- quan_vec(physeq)
otudb3<-make_quan(physeq,qp[1])
#the old filter
FSr  = transform_sample_counts(physeq, function(x) x / sum(x) )
otudb3 = filter_taxa(FSr, function(x) sum(x) > .05, TRUE)
#fill graph, order the categories
sample_data(otudb3)$category_site<-factor(sample_data(otudb3)$category_site,levels = c('Control_stomach','EGS_stomach','Control_ileum','EGS_ileum','Control_caecum','EGS_caecum','Control_colon','EGS_colon','Control_faeces','EGS_faeces','Co-Grazing_faeces'))
sample_data(otudb3)$site<-factor(sample_data(otudb3)$site,levels = c('stomach','ileum','caecum','colon','faeces'))
sample_data(physeq)$category_site<-factor(sample_data(physeq)$category_site,levels = c('Control_stomach','EGS_stomach','Control_ileum','EGS_ileum','Control_caecum','EGS_caecum','Control_colon','EGS_colon','Control_faeces','EGS_faeces','Co-Grazing_faeces'))
sample_data(physeq)$site<-factor(sample_data(physeq)$site,levels = c('stomach','ileum','caecum','colon','faeces'))
#color vector
mycolors<-c(brewer.pal(12, "Paired"),brewer.pal(12, "Set3"),brewer.pal(12, "Set3"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"),brewer.pal(8, "Pastel2"),brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"),brewer.pal(9, "Pastel1"), wes_palettes$"BottleRocket1",wes_palettes$"BottleRocket2",wes_palettes$"Rushmore1",wes_palettes$"Rushmore",wes_palettes$"Royal1",wes_palettes$"Royal2",wes_palettes$"Zissou1",wes_palettes$"Darjeeling1",wes_palettes$"Darjeeling2",wes_palettes$"Chevalier1",wes_palettes$"FantasticFox1",wes_palettes$"Moonrise1",wes_palettes$"Moonrise2",wes_palettes$"Moonrise3",wes_palettes$"Cavalcanti1",wes_palettes$"GrandBudapest1",wes_palettes$"GrandBudapest2",wes_palettes$"IsleofDogs1",wes_palettes$"IsleofDogs2")
#a function for fill-kingdom&phylum
plot_bar2 <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL)
{
    mdf = psmelt(physeq)
    p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
    p = p + geom_bar(stat = "identity", position = "fill", color = "black")#change position to 'fill'
    p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
    if (!is.null(facet_grid)) {
        p <- p + facet_grid(facet_grid)
    }
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    return(p)
}

#make & save the plot, unfiltered, seperate the egs and the control,color is managed
#kingdom
#filtered
egsvsctrl_kf <- plot_bar2(tax_glom(otudb3, taxrank=rank_names(otudb3)[2], NArm=TRUE),x='category_site',fill='K',title='Abundance by Kingdom filtered(the total read count of an OTU > 0.05 of all OTUs)') +geom_bar(aes(x=category_site,color=K, fill=K), stat='identity', position="fill")+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
ggsave(filename='egsvsctrl_k_filtered.pdf',plot=egsvsctrl_kf,device='pdf',dpi=300,width=8)
#unfiltered
egsvsctrl_k <- plot_bar2(tax_glom(physeq, taxrank=rank_names(physeq)[2], NArm=TRUE),x='category_site',fill='K',title='Abundance by Kingdom') +geom_bar(aes(x=category_site,color=K, fill=K), stat='identity', position="fill")+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
ggsave(filename='egsvsctrl_k.pdf',plot=egsvsctrl_k,device='pdf',dpi=300)
#phylum
#filtered
egsvsctrl_pf <- plot_bar2(tax_glom(otudb3, taxrank=rank_names(otudb3)[2], NArm=TRUE),x='category_site',fill='P',title='Abundance by Phylum filtered(the total read count of an OTU > 0.05 of all OTUs)') +geom_bar(aes(x=category_site,color=P, fill=P), stat='identity', position="fill")+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
ggsave(filename='egsvsctrl_p_filtered.pdf',plot=egsvsctrl_pf,device='pdf',dpi=300,width=8)
#unfiltered
egsvsctrl_p <- plot_bar2(tax_glom(physeq, taxrank=rank_names(physeq)[2], NArm=TRUE),x='category_site',fill='P',title='Abundance by Phylum') +geom_bar(aes(x=category_site,color=P, fill=P), stat='identity', position="fill")+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
ggsave(filename='egsvsctrl_p.pdf',plot=egsvsctrl_p,device='pdf',dpi=300)

#other graphs needs reordering the legends:
#SOME WORK IN LINUX
sort -t$',' -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 taxa.txt|cut -d ',' -f 4|uniq -> class.txt
sort -t$',' -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 taxa.txt|cut -d ',' -f 5|uniq -> order.txt
sort -t$',' -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 taxa.txt|cut -d ',' -f 6|uniq -> family.txt
#remove extra 'unidentified' (add back at the end manually)
sed '/^unidentified/d' class.txt |sed '/^NA/d'>C.txt
sed '/^unidentified/d' order.txt |sed '/^NA/d'>O.txt
sed '/^unidentified/d' family.txt |sed '/^NA/d'>F.txt
#use nano to 1: remove C,O,F, 2, add unidentified, 3,replicate the first line!!!
#another way is to use c(read.csv('C.txt'))$C)and put 'C' as the first line
#for class:
fill_graph<-function(){
#class filtered
x<-tax_glom(otudb3, taxrank=rank_names(otudb3)[3])
xmelt<-psmelt(x)
xmelt$C<-factor(xmelt$C,levels=c(read.csv('C.txt'))$C)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='C'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=c(mycolors,mycolors))+ scale_colour_manual(values=c(mycolors,mycolors))
xp=xp+ggtitle('Abundance by Class from OTUs filtered(the total read count of an OTU > 0.05 of all OTUs)')
xp=xp+geom_bar(aes(x=category_site,color=C, fill=C), stat='identity', position="fill")
ggsave(filename='egsvsctrl_c_filtered.pdf',plot=xp,device='pdf',width=8,height=6,dpi=300)
#class unfiltered
x<-tax_glom(physeq, taxrank=rank_names(physeq)[3])
xmelt<-psmelt(x)
xmelt$C<-factor(xmelt$C,levels=c(read.csv('C.txt'))$C)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='C'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=c(mycolors,mycolors))+ scale_colour_manual(values=c(mycolors,mycolors))
xp=xp+ggtitle('Abundance by Class from OTUs, unfiltered')
xp=xp+geom_bar(aes(x=category_site,color=C, fill=C), stat='identity', position="fill")
ggsave(filename='egsvsctrl_c.pdf',plot=xp,device='pdf',width=20,height=6,dpi=300)
#order filtered
x<-tax_glom(otudb3, taxrank=rank_names(otudb3)[4])
xmelt<-psmelt(x)
xmelt$O<-factor(xmelt$O,levels=c(read.csv('O.txt'))$O)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='O'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=c(mycolors,mycolors))+ scale_colour_manual(values=c(mycolors,mycolors))
xp=xp+ggtitle('Abundance by Order from OTUs filtered(the total read count of an OTU > 0.05 of all OTUs)')
xp=xp+geom_bar(aes(x=category_site,color=O, fill=O), stat='identity', position="fill")
ggsave(filename='egsvsctrl_o_filtered.pdf',plot=xp,device='pdf',width=10,height=6,dpi=300)
#order unfiltered
x<-tax_glom(physeq, taxrank=rank_names(physeq)[4])
xmelt<-psmelt(x)
xmelt$O<-factor(xmelt$O,levels=c(read.csv('O.txt'))$O)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='O'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=c(mycolors,mycolors))+ scale_colour_manual(values=c(mycolors,mycolors))
xp=xp+ggtitle('Abundance by Order from OTUs, unfiltered')
xp=xp+geom_bar(aes(x=category_site,color=O, fill=O), stat='identity', position="fill")
ggsave(filename='egsvsctrl_o.pdf',plot=xp,device='pdf',width=25,height=6,dpi=300)
#family filtered
x<-tax_glom(otudb3, taxrank=rank_names(otudb3)[5])
xmelt<-psmelt(x)
xmelt$F<-factor(xmelt$F,levels=c(read.csv('F.txt'))$F)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='F'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=c(mycolors,mycolors))+ scale_colour_manual(values=c(mycolors,mycolors))
xp=xp+ggtitle('Abundance by family from OTUs filtered(the total read count of an OTU > 0.05 of all OTUs)')
xp=xp+geom_bar(aes(x=category_site,color=F, fill=F), stat='identity', position="fill")
ggsave(filename='egsvsctrl_f_filtered.pdf',plot=xp,device='pdf',width=15,height=7,dpi=300)
#family unfiltered
x<-tax_glom(physeq, taxrank=rank_names(physeq)[5])
xmelt<-psmelt(x)
xmelt$F<-factor(xmelt$F,levels=c(read.csv('F.txt'))$F)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='F'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=c(mycolors,mycolors))+ scale_colour_manual(values=c(mycolors,mycolors))
xp=xp+ggtitle('Abundance by family from OTUs, unfiltered')
xp=xp+geom_bar(aes(x=category_site,color=F, fill=F), stat='identity', position="fill")
ggsave(filename='egsvsctrl_f.pdf',plot=xp,device='pdf',width=49,height=7,dpi=300)
}
fill_graph()
#
#
#
#
#
#
#
#
#now for stack
#define a new function plot_bar3,change position fill to stack
plot_bar3 <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL)
{
    mdf = psmelt(physeq)
    p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
    p = p + geom_bar(stat = "identity", position = "stack", color = "black")#change position to 'fill'
    p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
    if (!is.null(facet_grid)) {
        p <- p + facet_grid(facet_grid)
    }
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    return(p)
}

#make & save the plot, unfiltered, seperate the egs and the control,color is managed

#kingdom
#filtered
egsvsctrl_kf <- plot_bar3(tax_glom(otudb3, taxrank=rank_names(otudb3)[2], NArm=TRUE),x='category_site',fill='K',title='Abundance by Kingdom filtered(the total read count of an OTU > 0.05 of all OTUs)') +geom_bar(aes(x=category_site,color=K, fill=K), stat='identity', position="stack")+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
ggsave(filename='egsvsctrl_K_filtered.pdf',plot=egsvsctrl_kf,device='pdf',dpi=300,width=8)
#unfiltered
egsvsctrl_k <- plot_bar3(tax_glom(physeq, taxrank=rank_names(physeq)[2], NArm=TRUE),x='category_site',fill='K',title='Abundance by Kingdom') +geom_bar(aes(x=category_site,color=K, fill=K), stat='identity', position="stack")+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
ggsave(filename='egsvsctrl_k.pdf',plot=egsvsctrl_k,device='pdf',dpi=300)

#phylum
#filtered
egsvsctrl_pf <- plot_bar3(tax_glom(otudb3, taxrank=rank_names(otudb3)[2], NArm=TRUE),x='category_site',fill='P',title='Abundance by Phylum filtered (the sum of an OTU across all samples > 0.05 of all OTUs)') +geom_bar(aes(x=category_site,color=P, fill=P), stat='identity', position="stack")+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
ggsave(filename='egsvsctrl_p_filtered.pdf',plot=egsvsctrl_pf,device='pdf',dpi=300,width=8)
#unfiltered
egsvsctrl_p <- plot_bar3(tax_glom(physeq, taxrank=rank_names(physeq)[2], NArm=TRUE),x='category_site',fill='P',title='Abundance by Phylum') +geom_bar(aes(x=category_site,color=P, fill=P), stat='identity', position="stack")+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
ggsave(filename='egsvsctrl_p.pdf',plot=egsvsctrl_p,device='pdf',dpi=300)
#class filtered
x<-tax_glom(otudb3, taxrank=rank_names(otudb3)[3])
xmelt<-psmelt(x)
xmelt$C<-factor(xmelt$C,levels=c(read.csv('C.txt'))$C)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='C'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
xp=xp+ggtitle('Abundance by Class filtered (top 50% quantile)')
xp=xp+geom_bar(aes(x=category_site,color=C, fill=C), stat='identity', position="stack")
ggsave(filename='egsvsctrl_c_filtered.pdf',plot=xp,device='pdf',width=8,height=7,dpi=300)
#unfiltered
x<-tax_glom(physeq, taxrank=rank_names(physeq)[3])
xmelt<-psmelt(x)
xmelt$C<-factor(xmelt$C,levels=c(read.csv('C.txt'))$C)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='C'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
xp=xp+ggtitle('Abundance by Class from OTUs, unfiltered')
xp=xp+geom_bar(aes(x=category_site,color=C, fill=C), stat='identity', position="stack")
ggsave(filename='egsvsctrl_c.pdf',plot=xp,device='pdf',width=10,height=6,dpi=300)
#for order
x<-tax_glom(otudb3, taxrank=rank_names(otudb3)[4])
xmelt<-psmelt(x)
xmelt$O<-factor(xmelt$O,levels=c(read.csv('O.txt'))$O)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='O'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
xp=xp+ggtitle('Abundance by Order filtered (top 50% quantile)')
xp=xp+geom_bar(aes(x=category_site,color=O, fill=O), stat='identity', position="stack")
ggsave(filename='egsvsctrl_o_filtered.pdf',plot=xp,device='pdf',width=10,height=6,dpi=300)
#unfiltered
x<-tax_glom(physeq, taxrank=rank_names(physeq)[4])
xmelt<-psmelt(x)
xmelt$O<-factor(xmelt$O,levels=c(read.csv('O.txt'))$O)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='O'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
xp=xp+ggtitle('Abundance by Order from OTUs, unfiltered')
xp=xp+geom_bar(aes(x=category_site,color=O, fill=O), stat='identity', position="stack")
ggsave(filename='egsvsctrl_o.pdf',plot=xp,device='pdf',width=25,height=7,dpi=300)
#for family:
x<-tax_glom(otudb3, taxrank=rank_names(otudb3)[5])
xmelt<-psmelt(x)
xmelt$F<-factor(xmelt$F,levels=c(read.csv('F.txt'))$F)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='F'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=mycolors)+ scale_colour_manual(values=mycolors)
xp=xp+ggtitle('Abundance by family from OTUs, filtered (top 50% quantile)')
xp=xp+geom_bar(aes(x=category_site,color=F, fill=F), stat='identity', position="stack")
ggsave(filename='egsvsctrl_f_filtered.pdf',plot=xp,device='pdf',width=15,height=6,dpi=300)
#unfiltered
x<-tax_glom(physeq, taxrank=rank_names(physeq)[5])
xmelt<-psmelt(x)
xmelt$F<-factor(xmelt$F,levels=c(read.csv('F.txt'))$F)
xp=ggplot(xmelt,aes_string(x='category_site',y='Abundance',fill='F'))+ theme(axis.text.x = element_text(angle = 90))
xp=xp+ scale_fill_manual(values=c(mycolors,mycolors))+ scale_colour_manual(values=c(mycolors,mycolors))
xp=xp+ggtitle('Abundance by family from OTUs, unfiltered')
xp=xp+geom_bar(aes(x=category_site,color=F, fill=F), stat='identity', position="stack")
ggsave(filename='egsvsctrl_f.pdf',plot=xp,device='pdf',width=25,height=6,dpi=300)





















