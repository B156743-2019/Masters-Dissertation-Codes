call_library <- function() {
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
library('mixOmics')
library("wesanderson")
library("writexl")
library('pheatmap')}
#make tax_table
make_tax_table<- function(inputfile){ 
taxa_otu <- read.csv(inputfile)
taxa_otu<-unique(taxa_otu)
rownames(taxa_otu) <- taxa_otu[1:nrow(taxa_otu),1:1]
taxa_otu <- taxa_otu[1:nrow(taxa_otu),2:8]
taxa_otu$K<-gsub("ssignable","unidentified",taxa_otu$K)
taxa_otu$K<-gsub("_Fungi","Fungi",taxa_otu$K)
taxa_otu<-as.matrix(taxa_otu)
taxa_otu<-tax_table(taxa_otu)
return (taxa_otu)
}
taxa_otu<-make_tax_table('db2otutaxa.txt')
taxa<-make_tax_table('taxa.txt')
#intersections
#in new but not in old:
taxa_set <- function (set_operation,x,y,filename){
inter<-x[set_operation(rownames(x),rownames(y))]
df<-data.frame(OTU=rownames(inter),K=inter[,1],P=inter[,2],C=inter[,3],O=inter[,4],F=inter[,5],G=inter[,6],S=inter[,7])
#write_xlsx(df,filename)
return(df)
}
#excution
inter<-taxa_set(intersect,taxa,taxa_otu,'intersect_old_new.xls')
#taxa_set(setdiff,taxa,taxa_otu,'setdiff_new_old.xls')
diff<-taxa_set(setdiff,taxa_otu,taxa,'setdiff_old_new.xls')
#to merge the two sets together
taxa_otu<-data.frame(
OTU=append(as.character(inter$OTU),as.character(diff$OTU)),K=append(as.character(inter$K),as.character(diff$K)),
P=append(as.character(inter$P),as.character(diff$P)),C=append(as.character(inter$C),as.character(diff$C)),
O=append(as.character(inter$O),as.character(diff$O)),F=append(as.character(inter$F),as.character(diff$F)),
G=append(as.character(inter$G),as.character(diff$G)),S=append(as.character(inter$S),as.character(diff$S))
)
#remake the tax_file
rownames(taxa_otu) <- taxa_otu[1:nrow(taxa_otu),1:1]
taxa_otu <- taxa_otu[1:nrow(taxa_otu),2:8]
taxa_otu$K<-gsub("ssignable","unidentified",taxa_otu$K)
taxa_otu$K<-gsub("_Fungi","Fungi",taxa_otu$K)
taxa_otu<-as.matrix(taxa_otu)
taxa_otu<-tax_table(taxa_otu)
#copy the physeq1 data from previous analysis, change taxa table, and collapse by species
physeq1 = phyloseq(otu_table(physeq1), taxa_otu, sample_data(physeq1), phy_tree(physeq1))
by_K<-tax_glom(physeq1, taxrank=rank_names(physeq1)[1])
by_P<-tax_glom(physeq1, taxrank=rank_names(physeq1)[2])
by_C<-tax_glom(physeq1, taxrank=rank_names(physeq1)[3])
by_O<-tax_glom(physeq1, taxrank=rank_names(physeq1)[4])
by_F<-tax_glom(physeq1, taxrank=rank_names(physeq1)[5])
by_G<-tax_glom(physeq1, taxrank=rank_names(physeq1)[6])
by_S<-tax_glom(physeq1, taxrank=rank_names(physeq1)[7])
#the basic functions
#filter by prevalence
site_prevalence<- function(physeqobj,site_string) {
x<-otu_table(physeqobj)[,sample_data(physeqobj)$category_site==paste('EGS',site_string,sep='_')]
zeros<-data.frame(apply(x,1,function(y) length(which(y == 0))))
names(zeros)[1]<-'zcount'
prev<-phyloseq(
otu_table(physeqobj)[rownames(subset(zeros,zcount<ncol(x)/2))],
tax_table(physeqobj)[rownames(subset(zeros,zcount<ncol(x)/2))],
sample_data(physeqobj))
return (prev)
}
#filter by 50% quantile
quan_vec <- function (physeqobj) {
return (quantile(rowSums(otu_table(physeqobj)), probs = c(0.5), na.rm = FALSE, names = TRUE, type = 7))
}
#select the desired OTUs (after seeing the resulted quantiles):
make_quan <- function (physeqobj,cut_int){
quan<-phyloseq(
otu_table(physeqobj)[rowSums(otu_table(physeqobj))>cut_int],
tax_table(physeqobj)[rowSums(otu_table(physeqobj))>cut_int],
sample_data(physeqobj))
return (quan)
}
#DESeq2 analysis on one site
summary_table <- function (physeqobj,result,site_string,filename,category_str1,category_str2,level_int){
p_site<- lfcShrink(dds=result, contrast = c("category_site", paste(category_str1,site_string,sep='_'), paste(category_str2,site_string,sep='_')), type='ashr')
p_site<-subset(p_site,padj<0.05)
p_site<-data.frame(OTU=rownames(p_site),p_site,unlist(tax_table(physeqobj)[row.names(p_site),1:level_int]),site=c(site_string))
write_xlsx(p_site,filename)
}
#put the functions together
pre_quan_seq <- function (physeqobj, site_string,category_str1,category_str2,level=NULL,level_int){
#prev<-site_prevalence(physeqobj,site_string)
#qp <- quan_vec(physeqobj)
#quan<-make_quan(prev,qp[1])
quan<-physeqobj
otu_table(quan)[otu_table(quan)==0]<-1 #get rid of 1
dseqcs <-phyloseq_to_deseq2(quan, design = ~ category_site)
resultcs <- DESeq(dseqcs)
summary_table(physeqobj,resultcs,site_string,paste(level,category_str1,category_str2,site_string,'.xls',sep='_'),category_str1,category_str2,level_int)
}
#all sites
all_type <- function(physeqobj,category_str1,category_str2,level,level_int){
pre_quan_seq(physeqobj,'stomach',category_str1,category_str2,level,level_int)
pre_quan_seq(physeqobj,'ileum',category_str1,category_str2,level,level_int)
pre_quan_seq(physeqobj,'caecum',category_str1,category_str2,level,level_int)
pre_quan_seq(physeqobj,'colon',category_str1,category_str2,level,level_int)
pre_quan_seq(physeqobj,'faeces',category_str1,category_str2,level,level_int)
pre_quan_seq(physeqobj,'faeces','EGS','Co-Grazing',level,level_int)
pre_quan_seq(physeqobj,'faeces','Co-Grazing','Control',level,level_int)
}

#excution
all_type(physeq1,'EGS', 'Control', level='OTU',7)
#excution at each level
all_type(by_K,'EGS', 'Control',level='K',1)
all_type(by_P,'EGS', 'Control',level='P',2)
all_type(by_C,'EGS', 'Control',level='C',3)
all_type(by_O,'EGS', 'Control',level='O',4)
all_type(by_F,'EGS', 'Control',level='F',5)
all_type(by_G,'EGS', 'Control',level='G',6)
all_type(by_S,'EGS', 'Control',level='S',7)#done
#get heatmaps for EGS and Ctrl
heat_map <- function (physeqobj,result,site_string){
deseqsig<- lfcShrink(dds=result, contrast = c("category_site", paste('EGS',site_string,sep='_'), paste('Control',site_string,sep='_')), type='ashr')
deseqsig<-subset(deseqsig,padj<0.05)
deseqsig<-subset(deseqsig,log2FoldChange<0)
#start to add the pheatmap part
#db<-otu_table(physeqobj)[rownames(deseqsig),][,sample_data(physeqobj)$site==site_string]
x<-subset(sample_data(physeq1),category=='EGS'|category=='Control')
x<-subset(x,site==site_string)
db<-otu_table(physeqobj)[rownames(deseqsig),][,rownames(x)]
#db<-db[,colnames(db) != '1000']
db[db==0]<-0.00001
db<-apply(db,1, function(x) log10(x/sum(x)))
anno<-data.frame(category=sample_data(physeqobj)[sample_data(physeqobj)$site==site_string,]$category)
rownames(anno)<-rownames(sample_data(physeqobj)[sample_data(physeqobj)$site==site_string,])
pheatmap(t(db),annotation_col=anno,fontsize_row = 6,cellheight=5,main=paste(site_string,':','significantly underrepresented taxa in EGS horses', sep=' '))
}
#put the functions together
pre_quan_heapmap <- function (physeqobj, site_string){
prev<-site_prevalence(physeqobj,site_string)
#prev<-physeqobj
qp <- quan_vec(physeqobj)
quan<-make_quan(prev,qp[1])
otu_table(quan)[otu_table(quan)==0]<-1 #get rid of 1
dseqcs <-phyloseq_to_deseq2(quan, design = ~ category_site)
resultcs <- DESeq(dseqcs)
heat_map(physeqobj,resultcs,site_string)
}
#
all_heatmap <- function(physeqobj){
pdf('EGS_vs_Ctrl.pdf',height=40)
pre_quan_heapmap(physeqobj,'stomach')
pre_quan_heapmap(physeqobj,'ileum')
pre_quan_heapmap(physeqobj,'caecum')
pre_quan_heapmap(physeqobj,'colon')
pre_quan_heapmap(physeqobj,'faeces')
dev.off()
}
#excution
all_heatmap(physeq1)
#collapse to get the phylotype version of analysis
all_heatmap(by_S)










#comparison by categories, OTU and Phylotypes, unfiltered
#category prevalence
cate_prevalence<- function(physeqobj) {
x<-otu_table(physeqobj)[,sample_data(physeqobj)$category=='EGS']
zeros<-data.frame(apply(x,1,function(y) length(which(y == 0))))
names(zeros)[1]<-'zcount'
prev<-phyloseq(
otu_table(physeqobj)[rownames(subset(zeros,zcount<ncol(x)/2))],#that is prevalence>50%, or less than 50% zero count
tax_table(physeqobj)[rownames(subset(zeros,zcount<ncol(x)/2))],
sample_data(physeqobj))
return (prev)
}
#
summary_table_cat <- function (physeqobj,result,category_str1,category_str2,type_str,type_int){
p_site<- lfcShrink(dds=result, contrast = c("category", category_str1, category_str2), type='ashr')
p_site<-subset(p_site,padj<0.05)
p_site<-data.frame(OTU=rownames(p_site),p_site,unlist(tax_table(physeqobj)[row.names(p_site),1:type_int]))
write_xlsx(p_site,paste(type_str,category_str1,category_str2,'.xlsx',sep='_'))
}
#put the functions together
pre_quan_seq_cat <- function (physeqobj,type_str,type_int){
#prev<-cate_prevalence(physeqobj)
#qp <- quan_vec(physeqobj)
#quan<-make_quan(prev,qp[1])
quan<-physeqobj
otu_table(quan)[otu_table(quan)==0]<-1 #get rid of 1
dseqcs <-phyloseq_to_deseq2(quan, design = ~ category)
resultcs <- DESeq(dseqcs)
summary_table_cat(physeqobj,resultcs,'EGS','Control',type_str,type_int)
summary_table_cat(physeqobj,resultcs,'EGS','Co-Grazing',type_str,type_int)
summary_table_cat(physeqobj,resultcs,'Co-Grazing','Control',type_str,type_int)
}

#excution
pre_quan_seq_cat(physeq1,'OTU')
pre_quan_seq_cat(by_S,'phylotype')
pre_quan_seq_cat(by_P,'P',2)
pre_quan_seq_cat(by_C,'C',3)
pre_quan_seq_cat(by_O,'O',4)
pre_quan_seq_cat(by_F,'F',5)
pre_quan_seq_cat(by_G,'G',6)
































