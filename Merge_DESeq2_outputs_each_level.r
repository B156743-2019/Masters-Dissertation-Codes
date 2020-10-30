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
library('pheatmap')
library('readxl')}

load('0706.rdata')
#rename the files
rename_files<-function(level_str){
files <- list.files(pattern="*.xls")
newfiles <- gsub(".xls$", ".xlsx", files)
newfiles2<-gsub(level_str,'',newfiles)
file.rename(files, newfiles2)}

build_df<-function(rename_from,level_int,filename){
rename_files(rename_from)
gcf<-data.frame(read_excel("Co-Grazing_Control_faeces_.xlsx",sheet=1))
egf<-data.frame(read_excel("EGS_Co-Grazing_faeces_.xlsx",sheet=1))
ecca<-data.frame(read_excel("EGS_Control_caecum_.xlsx",sheet=1))
ecco<-data.frame(read_excel("EGS_Control_colon_.xlsx" ,sheet=1))
eci<-data.frame(read_excel("EGS_Control_ileum_.xlsx" ,sheet=1))
ecs<-data.frame(read_excel("EGS_Control_stomach_.xlsx",sheet=1))
ecf<-data.frame(read_excel("EGS_Control_faeces_.xlsx",sheet=1))
eco<-data.frame(read_excel("EGS_Control_.xlsx",sheet=1))

otus<-unique(union(union(union(union(union(union(union(gcf$OTU,egf$OTU),ecca$OTU),ecco$OTU),eci$OTU),ecs$OTU),eco$OTU),ecf$OTU))
merge_table<-data.frame(OTU=rownames(taxa_otu[otus,]),taxa_otu[otus,1:level_int])
#populating the dataframe
merge_table$egs_ctrl_ov_LFC<-eco$log2FoldChange[match(rownames(merge_table),eco$OTU)]
merge_table$egs_ctrl_ov_padj<-eco$padj[match(rownames(merge_table),eco$OTU)]

merge_table$egs_ctrl_st_LFC<-ecs$log2FoldChange[match(rownames(merge_table),ecs$OTU)]
merge_table$egs_ctrl_st_padj<-ecs$padj[match(rownames(merge_table),ecs$OTU)]

merge_table$egs_ctrl_il_LFC<-eci$log2FoldChange[match(rownames(merge_table),eci$OTU)]
merge_table$egs_ctrl_il_padj<-eci$padj[match(rownames(merge_table),eci$OTU)]

merge_table$egs_ctrl_ca_LFC<-ecca$log2FoldChange[match(rownames(merge_table),ecca$OTU)]
merge_table$egs_ctrl_ca_padj<-ecca$padj[match(rownames(merge_table),ecca$OTU)]

merge_table$egs_ctrl_co_LFC<-ecco$log2FoldChange[match(rownames(merge_table),ecco$OTU)]
merge_table$egs_ctrl_co_padj<-ecco$padj[match(rownames(merge_table),ecco$OTU)]

merge_table$egs_ctrl_fa_LFC<-ecf$log2FoldChange[match(rownames(merge_table),ecf$OTU)]
merge_table$egs_ctrl_fa_padj<-ecf$padj[match(rownames(merge_table),ecf$OTU)]

merge_table$egs_cog_fa_LFC<-egf$log2FoldChange[match(rownames(merge_table),egf$OTU)]
merge_table$egs_cog_fa_padj<-egf$padj[match(rownames(merge_table),egf$OTU)]

merge_table$cog_ctrl_fa_LFC<-gcf$log2FoldChange[match(rownames(merge_table),gcf$OTU)]
merge_table$cog_ctrl_fa_padj<-gcf$padj[match(rownames(merge_table),gcf$OTU)]

write_xlsx(merge_table,filename)
#return (merge_table)
}
build_df('*P_',2,'merged_phylum.xlsx')
build_df('*C_',3,'merged_class.xlsx')
build_df('*O_',4,'merged_order.xlsx')
build_df('*F_',5,'merged_family.xlsx')
build_df('*G_',6,'merged_genus.xlsx')

#use a loop instead
loop_build_df<-function(level_str,level_int,filename){
file_rename(level_str)
file_list<-c("EGS_Control_.xlsx","EGS_Control_stomach_.xlsx","EGS_Control_ileum_.xlsx",
"EGS_Control_caecum_.xlsx","EGS_Control_colon_.xlsx","EGS_Control_faeces_.xlsx",
"EGS_Co-Grazing_faeces_.xlsx")
up_to_now<-c()
for (i in 1:length(file_list)){
df<-data.frame(read_excel(file_list[i],sheet=1))
union_otu<-union(df$OTU,up_to_now)
up_to_now<-union_otu
}
otus<-up_to_now
merge_table<-data.frame(OTU=rownames(taxa_otu[otus,]),taxa_otu[otus,1:level_int])
name_list<- c("EGS_Control_overall","EGS_Control_stomach","EGS_Control_ileum",
"EGS_Control_caecum","EGS_Control_colon","EGS_Control_faeces","EGS_Co-Grazing_faeces")
for (i in 1:length(name_list)){
df<-data.frame(read_excel(file_list[i],sheet=1))
merge_table$LFC<-df$log2FoldChange[match(rownames(merge_table),df$OTU)]#to create the columns
names(merge_table)[names(merge_table) == "LFC"] <- paste(name_list[i],'LFC')#to rename the coloums
merge_table$padj<-df$padj[match(rownames(merge_table),df$OTU)]#to create the columns
names(merge_table)[names(merge_table) == "padj"] <-paste(name_list[i],'padj')#to rename the coloums
}
write_xlsx(merge_table,filename)
}

loop_build_df('*P_',2,'merged_phylum.xlsx')
loop_build_df('*C_',3,'merged_class.xlsx')
loop_build_df('*O_',4,'merged_order.xlsx')
loop_build_df('*F_',5,'merged_family.xlsx')
loop_build_df('*G_',6,'merged_genus.xlsx')







