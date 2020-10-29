#library
call_library<- function () {
library('scatterplot3d')
library('GUniFrac')
library('vegan')
library('ggplot2')
library('corrplot')
library('phyloseq')
library('devtools')
library('metagenomeSeq')
library('ape')
library('dplyr')
library('DESeq2') 
library('RVAideMemoire')
library("wesanderson")
library('mixOmics')
library('RVAideMemoire')
library("writexl")
}
#PLSDA
make_matrix<-function(physeqobj,site_string){
x<-otu_table(physeqobj)[,sample_data(physeqobj)$site==site_string]
x<-x[,colnames(x)!='1000']#remove fake sample
x<-subset(x,rowSums(x)!=0)
return(x)
}
#
make_factor<-function(physeqobj, site_string){
y<-sample_data(physeqobj)$category[sample_data(physeqobj)$site==site_string]
y<-subset(y,y!='Co-grazing')#remove fake sample
return(y)
}
#
make_names<-function(physeqobj,site_string) {
names<-rownames(sample_data(physeqobj)[sample_data(physeqobj)$site==site_string])
names<-subset(names,names!='1000')
return(names)
}
#
plsda_plot<-function(x,y,names,site_string,type_str){
plsda_xy<-splsda(t(x),y)
plotIndiv(plsda_xy, ind.names = NULL, ellipse = TRUE, legend =TRUE,title=paste('PLSDA Plot',site_string,'by',type_str,sep=' '))
#return(plsda_xy)
}
#
plsda_one_site<-function(physeqobj,site_string,type_str) {
x<-make_matrix(physeqobj,site_string)
y<-make_factor(physeqobj,site_string)
#names<-make_names(physeqobj,site_string)
plsda_plot(x,y,names,site_string,type_str)
}
#
plsda_all_sites<-function(physeqobj,file_name_string,type_str){
pdf(file_name_string,height=8)
plsda_one_site(physeqobj,'stomach',type_str)
plsda_one_site(physeqobj,'ileum',type_str)
plsda_one_site(physeqobj,'caecum',type_str)
plsda_one_site(physeqobj,'colon',type_str)
plsda_one_site(physeqobj,'faeces',type_str)
dev.off()
}
#excution
plsda_all_sites(physeq1,'PLSDA_OTUs.pdf','OTU')
plsda_all_sites(by_S,'PLSDA_phylotype.pdf','Phylotype')
#to show site differences in EGS samples
site_comparison<- function (physeqobj,category_string,type_string) {
x<-otu_table(physeqobj)[,sample_data(physeqobj)$category==category_string]
x<-x[,colnames(x)!='1000']#remove fake sample
x<-subset(x,rowSums(x)!=0)
y<-sample_data(physeqobj)$site[sample_data(physeqobj)$category==category_string]
y<-subset(y,y!='Co-grazing')#remove fake sample
names<-rownames(sample_data(physeqobj))[sample_data(physeqobj)$category==category_string]
names<-subset(names,names!='1000')
plsda_xy<-plsda(t(x),y)
plotIndiv(plsda_xy, ind.names = NULL, ellipse = TRUE, legend =TRUE,title=paste('Site Comparison of', category_string,'samples by',type_string,sep=' '))
}
#make the pdf
pdf('site_comparison.pdf',height=8,width=10)
site_comparison(physeq1,'EGS','OTU')
site_comparison(by_S,'EGS','phylotype')
site_comparison(physeq1,'Control','OTU')
site_comparison(by_S,'Control','phylotype')
dev.off()
#to show the difference between EGS and Control samples
category_comparison<- function (physeqobj,type_string,names_bolean=FALSE, cat_str1, cat_str2) {
x<-otu_table(physeqobj)[,(sample_data(physeqobj)$category==cat_str1)|(sample_data(physeqobj)$category==cat_str2)]
x<-x[,colnames(x)!='1000']#remove fake sample
x<-subset(x,rowSums(x)!=0)
y<-sample_data(physeqobj)$category[(sample_data(physeqobj)$category==cat_str1)|(sample_data(physeqobj)$category==cat_str2)]
y<-subset(y,y!='Co-grazing')#remove fake sample
names<-rownames(sample_data(physeqobj))[(sample_data(physeqobj)$category==cat_str1)|(sample_data(physeqobj)$category==cat_str2)]
names<-subset(names,names!='1000')
if(names_bolean=='TRUE'){names=names}#Show sample names on the graph
else{names='NULL'}
plsda_xy<-plsda(t(x),y)
plotIndiv(plsda_xy, ind.names = names, ellipse = TRUE, legend =TRUE,title=paste('Comparison of',cat_str1,'and',cat_str2,'samples by',type_string,sep=' '))
}
#make the pdf
pdf('Category_comparison.pdf',height=8,width=10)
category_comparison(physeq1,'OTU',names_bolean=TRUE,'EGS','Control')
category_comparison(by_S,'phylotype',names_bolean=TRUE,'EGS','Control')
category_comparison(physeq1,'OTU',names_bolean=FALSE,'EGS','Control')
category_comparison(by_S,'phylotype',names_bolean=FALSE,'EGS','Control')
category_comparison(physeq1,'OTU',names_bolean=TRUE,'EGS','Co-Grazing')
category_comparison(by_S,'phylotype',names_bolean=TRUE,'EGS','Co-Grazing')
category_comparison(physeq1,'OTU',names_bolean=FALSE,'EGS','Co-Grazing')
category_comparison(by_S,'phylotype',names_bolean=FALSE,'EGS','Co-Grazing')
category_comparison(physeq1,'OTU',names_bolean=TRUE,'Co-Grazing','Control')
category_comparison(by_S,'phylotype',names_bolean=TRUE,'Co-Grazing','Control')
category_comparison(physeq1,'OTU',names_bolean=FALSE,'Co-Grazing','Control')
category_comparison(by_S,'phylotype',names_bolean=FALSE,'Co-Grazing','Control')
dev.off()


#for all sites
plsda_vip<-function(physeqobj,file_name_string){
df<-data.frame(
vip_score(physeqobj,'stomach'),
vip_score(physeqobj,'ileum'),
vip_score(physeqobj,'caecum'),
vip_score(physeqobj,'colon'),
vip_score(physeqobj,'faeces'))
write_xlsx(df,file_name_string)
}
#excution
plsda_vip(physeq1,'vip_score_top_100_otu.xls')
plsda_vip(by_S,'vip_score_top_100_phylotype.xls')

#VIP score for overall comparison
vip_overall<- function (physeqobj,file_name_string) {
x<-otu_table(physeqobj)[,sample_data(physeqobj)$category!='Co-Grazing']
x<-x[,colnames(x)!='1000']#remove fake sample
x<-subset(x,rowSums(x)!=0)
y<-sample_data(physeqobj)$category[sample_data(physeqobj)$category!='Co-Grazing']
y<-subset(y,y!='Co-grazing')#remove fake sample
plsda_xy<-plsda(t(x),y)
df<-data.frame(
otu=rownames(head(PLSDA.VIP(plsda_xy)$tab,100)),
kingdom=unlist(tax_table(physeqobj)[rownames(head(PLSDA.VIP(plsda_xy)$tab,100))])[,1],
phylum=unlist(tax_table(physeqobj)[rownames(head(PLSDA.VIP(plsda_xy)$tab,100))])[,2],
class=unlist(tax_table(physeqobj)[rownames(head(PLSDA.VIP(plsda_xy)$tab,100))])[,3],
order=unlist(tax_table(physeqobj)[rownames(head(PLSDA.VIP(plsda_xy)$tab,100))])[,4],
family=unlist(tax_table(physeqobj)[rownames(head(PLSDA.VIP(plsda_xy)$tab,100))])[,5],
genus=unlist(tax_table(physeqobj)[rownames(head(PLSDA.VIP(plsda_xy)$tab,100))])[,6],
species=unlist(tax_table(physeqobj)[rownames(head(PLSDA.VIP(plsda_xy)$tab,100))])[,7],
vip=head(PLSDA.VIP(plsda_xy)$tab,100))
write_xlsx(df,file_name_string)
}
#excution
vip_overall(physeq1,'EGS_vs_Control_vip_score_top_100_otu.xls')
vip_overall(by_S,'EGS_vs_Control_vip_score_top_100_phylotype.xls')
#match vip score to the deseq2 results
#VIP score for all otu
vip_score<-function(physeqobj,site_string){
x<-make_matrix(physeqobj,site_string)
y<-make_factor(physeqobj,site_string)
plsda_xy<-splsda(t(x),y)
df<-data.frame(
site=c(site_string),
otu=rownames(PLSDA.VIP(plsda_xy)$tab),
vip=PLSDA.VIP(plsda_xy)$tab)
return(df)
}
#wrap everything in a huge function
#upload the corresponding deseq2 files, and name it by the site, i.e. stomach.csv
vip_deseq<-function(physeqobj, site_string){
dseqobj<-read.csv(paste(site_string,'.csv',sep=''))
rownames(dseqobj)<-dseqobj$OTU
vip<-vip_score(physeqobj, site_string)
inter<-intersect(rownames(dseqobj),rownames(vip))
df<-data.frame(dseqobj, vip_score=vip[inter,]$VIP)
write_xlsx(df,paste('EGS_vs_Ctrl_',site_string,'.xls', sep=''))
}
#excution uncollapsed, otu
vip_deseq(physeq1,'stomach')
vip_deseq(physeq1,'ileum')
vip_deseq(physeq1,'caecum')
vip_deseq(physeq1,'colon')
vip_deseq(physeq1,'faeces')
#excution collapsed, sp
vip_deseq(by_S,'stomach')
vip_deseq(by_S,'ileum')
vip_deseq(by_S,'caecum')
vip_deseq(by_S,'colon')
vip_deseq(by_S,'faeces')
#assess the clustering through unifrac distance and adonis, change $category to fit need
select_group<-function(physeqobj,site_string){
phy<-phyloseq(otu_table(physeqobj)[,sample_data(physeqobj)$site==site_string & ((sample_data(physeqobj)$category=='Control')|(sample_data(physeqobj)$category=='Co-Grazing'))], 
tax_table(physeqobj))
coal_tree = rcoal(ntaxa(phy), rooted=TRUE, tip.label=taxa_names(phy))#build a tree
phy<-phyloseq(otu_table(phy),tax_table(phy), 
sample_data(physeqobj)[sample_data(physeqobj)$site==site_string], 
phy_tree(physeqobj))#rebuild the physeq obj
return(phy)
}
assess_plsda<-function(physeqobj,site_string){
phy<-select_group(physeqobj,site_string)
metadf <- data.frame(sample_data(phy))
unifrac.dist <- UniFrac(phy, weighted = TRUE,normalized = TRUE,parallel = FALSE,fast = TRUE)
permanova <- adonis(unifrac.dist ~ category, data = metadf)
print(permanova)
}
#clustering assessment: is it different enough: uncollapsed
assess_plsda(physeq1,'stomach')
assess_plsda(physeq1,'ileum')
assess_plsda(physeq1,'caecum')
assess_plsda(physeq1,'colon')
assess_plsda(physeq1,'faeces')
#clustering assessment: is it different enough: collapsed
assess_plsda(by_S,'stomach')
assess_plsda(by_S,'ileum')
assess_plsda(by_S,'caecum')
assess_plsda(by_S,'colon')
assess_plsda(by_S,'faeces')
#compare sites of each category
assess_cat_site<-function(physeqobj,site_string1,site_string2,cat_string){
phy<-phyloseq(
otu_table(physeqobj)[(sample_data(physeqobj)$site==site_string1|sample_data(physeqobj)$site==site_string2)& sample_data(physeqobj)$category==cat_string], 
tax_table(physeqobj))
coal_tree = rcoal(ntaxa(phy), rooted=TRUE, tip.label=taxa_names(phy))#build a tree
phy<-phyloseq(otu_table(phy),tax_table(phy), 
sample_data(physeqobj)[(sample_data(physeqobj)$site==site_string1|sample_data(physeqobj)$site==site_string2)& sample_data(physeqobj)$category==cat_string], 
phy_tree(physeqobj))
metadf <- data.frame(sample_data(phy))
unifrac.dist <- UniFrac(phy, weighted = TRUE,normalized = TRUE,parallel = FALSE,fast = TRUE)
permanova <- adonis(unifrac.dist ~ site, data = metadf)
print(permanova)
}
#EGS uncollapsed
assess_cat_site(physeq1,'stomach','ileum','EGS')
assess_cat_site(physeq1,'stomach','caecum','EGS')
assess_cat_site(physeq1,'stomach','colon','EGS')
assess_cat_site(physeq1,'stomach','faeces','EGS')
assess_cat_site(physeq1,'ileum','caecum','EGS')
assess_cat_site(physeq1,'ileum','colon','EGS')
assess_cat_site(physeq1,'ileum','faeces','EGS')
assess_cat_site(physeq1,'caecum','colon','EGS')
assess_cat_site(physeq1,'caecum','faeces','EGS')
assess_cat_site(physeq1,'faeces','colon','EGS')
#EGS collapsed
assess_cat_site(by_S,'stomach','ileum','EGS')
assess_cat_site(by_S,'stomach','caecum','EGS')
assess_cat_site(by_S,'stomach','colon','EGS')
assess_cat_site(by_S,'stomach','faeces','EGS')
assess_cat_site(by_S,'ileum','caecum','EGS')
assess_cat_site(by_S,'ileum','colon','EGS')
assess_cat_site(by_S,'ileum','faeces','EGS')
assess_cat_site(by_S,'caecum','colon','EGS')
assess_cat_site(by_S,'caecum','faeces','EGS')
assess_cat_site(by_S,'faeces','colon','EGS')
#Ctrl uncollapsed
assess_cat_site(physeq1,'stomach','ileum','Control')
assess_cat_site(physeq1,'stomach','caecum','Control')
assess_cat_site(physeq1,'stomach','colon','Control')
assess_cat_site(physeq1,'stomach','faeces','Control')
assess_cat_site(physeq1,'ileum','caecum','Control')
assess_cat_site(physeq1,'ileum','colon','Control')
assess_cat_site(physeq1,'ileum','faeces','Control')
assess_cat_site(physeq1,'caecum','colon','Control')
assess_cat_site(physeq1,'caecum','faeces','Control')
assess_cat_site(physeq1,'faeces','colon','Control')
#Ctrl collapsed
assess_cat_site(by_S,'stomach','ileum','Control')
assess_cat_site(by_S,'stomach','caecum','Control')
assess_cat_site(by_S,'stomach','colon','Control')
assess_cat_site(by_S,'stomach','faeces','Control')
assess_cat_site(by_S,'ileum','caecum','Control')
assess_cat_site(by_S,'ileum','colon','Control')
assess_cat_site(by_S,'ileum','faeces','Control')
assess_cat_site(by_S,'caecum','colon','Control')
assess_cat_site(by_S,'caecum','faeces','Control')
assess_cat_site(by_S,'faeces','colon','Control')
#overall EGS vs Control
select_cat <-function(physeqobj){
phy<-phyloseq(otu_table(physeqobj)[,(sample_data(physeqobj)$category=='Co-Grazing')|(sample_data(physeqobj)$category=='Control')], 
tax_table(physeqobj))
coal_tree = rcoal(ntaxa(phy), rooted=TRUE, tip.label=taxa_names(phy))#build a tree
phy<-phyloseq(otu_table(phy),tax_table(phy), 
sample_data(physeqobj)[(sample_data(physeqobj)$category=='Co-Grazing')|(sample_data(physeqobj)$category=='Control')], 
phy_tree(physeqobj))#rebuild the physeq obj
return(phy)
}
assess_cat<-function(physeqobj){
phy<-select_cat(physeqobj)
metadf <- data.frame(sample_data(phy))
unifrac.dist <- UniFrac(phy, weighted = TRUE,normalized = TRUE,parallel = FALSE,fast = TRUE)
permanova <- adonis(unifrac.dist ~ category, data = metadf)
print(permanova)
}
#excution
assess_cat(physeq1)
assess_cat(by_S)



