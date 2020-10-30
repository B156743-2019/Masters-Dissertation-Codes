#
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
#remove the fake sample
trim_fake<- function (physeqobj){
physeq<-phyloseq(
otu_table(physeqobj)[,colnames(otu_table(physeqobj))!='1000'],
sample_data(physeqobj)[rownames(sample_data(physeqobj))!='1000'],
tax_table(physeqobj),
phy_tree(physeqobj))
sample_data(physeq)$category_site<-factor(sample_data(physeq)$category_site,levels = c('Control_stomach','EGS_stomach','Control_ileum','EGS_ileum','Control_caecum','EGS_caecum','Control_colon','EGS_colon','Control_faeces','EGS_faeces','Co-Grazing_faeces'))
sample_data(physeq)$site<-factor(sample_data(physeq)$site,levels = c('stomach','ileum','caecum','colon','faeces'))
return(physeq)
}
#make accumulation plot
accm<- function () {
sp1 <- specaccum(t(otu_table(physeq1)))
#sp2 <- specaccum(t(otu_table(physeq1)), "random")
pdf('accumulation_plot2.pdf')
plot(sp1, ci=2, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
#boxplot(sp2, col="yellow", add=TRUE, pch="+")
dev.off()
}

#bar plot of monthly distribution
pdf('monthly_distribution.pdf')
barplot(month,xlab='month',col=c("darkblue","red"),names.arg=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),beside=TRUE,legend=c('EGS','Ctrl'))
dev.off()

#modify the extant plot_richness function to make box plot
alpha_box<-function(physeqobj,plot_type_string){
p<-plot_richness(trim_fake(physeqobj), measures=c(plot_type_string),x='site',color='category')
p$layers <- p$layers[-1]
p$layers <- p$layers[-1]#remove unwanted layers (the dots)
p+geom_boxplot() #use a boxplot instead
}
#make alpha diversity plot, dot plot and box plot
pdf('alpha_diversity.pdf',height=8,width=10)
plot_richness(trim_fake(physeq1), measures=c("Observed"),x='site',color='category')
alpha_box(physeq1,'Observed')
plot_richness(trim_fake(physeq1), measures=c("Chao1"),x='site',color='category')
alpha_box(physeq1,'Chao1')
plot_richness(trim_fake(physeq1), measures=c("ACE"),x='site',color='category')
alpha_box(physeq1,'ACE')
plot_richness(trim_fake(physeq1), measures=c("Shannon"),x='site',color='category')
alpha_box(physeq1,'Shannon')
plot_richness(trim_fake(physeq1), measures=c("Simpson"),x='site',color='category')
alpha_box(physeq1,'Simpson')
plot_richness(trim_fake(physeq1), measures=c("InvSimpson"),x='site',color='category')
alpha_box(physeq1,'InvSimpson')
plot_richness(trim_fake(physeq1), measures=c("Fisher"),x='site',color='category')
alpha_box(physeq1,'Fisher')
dev.off()
#statistics of alpha diversity
#Pr > F – This is the p-value associated with the F statistic of a given effect and test statistic. The null hypothesis that a given predictor has no effect on either of the outcomes is evaluated with regard to this p-value. For a given alpha level, if the p-value is less than alpha, the null hypothesis is rejected.
summary_cat<-function(){
print('#ANOVA between categories')
print(summary(aov(DF$Observed ~ DF$category)))
print(summary(aov(DF$Chao1 ~ DF$category)))
print(summary(aov(DF$ACE ~ DF$category)))
print(summary(aov(DF$Shannon ~ DF$category)))
print(summary(aov(DF$Simpson ~ DF$category)))
print(summary(aov(DF$InvSimpson ~ DF$category)))
print(summary(aov(DF$Fisher ~ DF$category)))
}
#summary_cat() #for comparison between categories using different parameters, all sites together
#for site
summary_site<-function(method){
alpha_obj<-0
if (method=='Observed') {alpha_obj<-DF$Observed}
else if (method=='Chao1') {alpha_obj<-DF$Chao1}
else if (method=='ACE') {alpha_obj<-DF$ACE}
else if (method=='Shannon') {alpha_obj<-DF$Shannon}
else if (method=='Simpson') {alpha_obj<-DF$Simpson}
else if(method=='InvSimpson') {alpha_obj<-DF$InvSimpson}
else if(method=='Fisher') {alpha_obj<-DF$Fisher}
print(paste('#ANOVA between sites using',method,sep=' '))
print(summary(aov(alpha_obj ~ DF$site)))
print(summary(aov(alpha_obj[(DF$site=='stomach')|(DF$site=='ileum')] ~ DF$site[(DF$site=='stomach')|(DF$site=='ileum')])))
print(summary(aov(alpha_obj[(DF$site=='stomach')|(DF$site=='caecum')] ~ DF$site[(DF$site=='stomach')|(DF$site=='caecum')])))
print(summary(aov(alpha_obj[(DF$site=='stomach')|(DF$site=='colon')] ~ DF$site[(DF$site=='stomach')|(DF$site=='colon')])))
print(summary(aov(alpha_obj[(DF$site=='stomach')|(DF$site=='faeces')] ~ DF$site[(DF$site=='stomach')|(DF$site=='faeces')])))
print(summary(aov(alpha_obj[(DF$site=='ileum')|(DF$site=='caecum')] ~ DF$site[(DF$site=='ileum')|(DF$site=='caecum')])))
print(summary(aov(alpha_obj[(DF$site=='ileum')|(DF$site=='colon')] ~ DF$site[(DF$site=='ileum')|(DF$site=='colon')])))
print(summary(aov(alpha_obj[(DF$site=='ileum')|(DF$site=='faeces')] ~ DF$site[(DF$site=='ileum')|(DF$site=='faeces')])))
print(summary(aov(alpha_obj[(DF$site=='caecum')|(DF$site=='colon')] ~ DF$site[(DF$site=='caecum')|(DF$site=='colon')])))
print(summary(aov(alpha_obj[(DF$site=='caecum')|(DF$site=='faeces')] ~ DF$site[(DF$site=='caecum')|(DF$site=='faeces')])))
print(summary(aov(alpha_obj[(DF$site=='faeces')|(DF$site=='colon')] ~ DF$site[(DF$site=='faeces')|(DF$site=='colon')])))
}
#category at each site
summary_cat_site<-function(method) {
alpha_obj<-0
if (method=='Observed') {alpha_obj<-DF$Observed}
else if (method=='Chao1') {alpha_obj<-DF$Chao1}
else if (method=='ACE') {alpha_obj<-DF$ACE}
else if (method=='Shannon') {alpha_obj<-DF$Shannon}
else if (method=='Simpson') {alpha_obj<-DF$Simpson}
else if(method=='InvSimpson') {alpha_obj<-DF$InvSimpson}
else if(method=='Fisher') {alpha_obj<-DF$Fisher}
print(paste('#ANOVA between EGS, Control, Co-Grazing at each site, using',method,sep=' '))
print(summary(aov(alpha_obj[DF$site=='stomach'] ~ DF$category[DF$site=='stomach'])))
print(summary(aov(alpha_obj[DF$site=='ileum'] ~ DF$category[DF$site=='ileum'])))
print(summary(aov(alpha_obj[DF$site=='caecum'] ~ DF$category[DF$site=='caecum'])))
print(summary(aov(alpha_obj[DF$site=='colon'] ~ DF$category[DF$site=='colon'])))
print(summary(aov(alpha_obj[DF$site=='faeces'] ~ DF$category[DF$site=='faeces'])))
}
#all together
alpha_stats<-function (physeqobj){
alpha<-data.frame(estimate_richness(physeqobj))
DF<- data.frame(sample_data(physeqobj),alpha)
summary_cat()
#compare sites
summary_site('Observed')
summary_site('Chao1')
summary_site('ACE')
summary_site('Shannon')
summary_site('Simpson')
summary_site('InvSimpson')
summary_site('Fisher')
#category at each site
summary_cat_site('Observed')
summary_cat_site('Chao1')
summary_cat_site('ACE')
summary_cat_site('Shannon')
summary_cat_site('Simpson')
summary_cat_site('InvSimpson')
summary_cat_site('Fisher')
}

#check who has higher alpha

mean(DF$Observed[(DF$site=='caecum')&(DF$category=='EGS')])/mean(DF$Observed[(DF$site=='caecum')&(DF$category=='Control')])
mean(DF$Observed[(DF$site=='colon')&(DF$category=='EGS')])/mean(DF$Observed[(DF$site=='colon')&(DF$category=='Control')])
mean(DF$Chao1[(DF$site=='caecum')&(DF$category=='EGS')])/mean(DF$Chao1[(DF$site=='caecum')&(DF$category=='Control')])
mean(DF$Chao1[(DF$site=='colon')&(DF$category=='EGS')])/mean(DF$Chao1[(DF$site=='colon')&(DF$category=='Control')])
mean(DF$ACE[(DF$site=='caecum')&(DF$category=='EGS')])/mean(DF$ACE[(DF$site=='caecum')&(DF$category=='Control')])
mean(DF$ACE[(DF$site=='colon')&(DF$category=='EGS')])/mean(DF$ACE[(DF$site=='colon')&(DF$category=='Control')])
mean(DF$Shannon[(DF$site=='faeces')&(DF$category=='EGS')])/mean(DF$Shannon[(DF$site=='faeces')&(DF$category=='Control')])
mean(DF$Shannon[(DF$site=='colon')&(DF$category=='EGS')])/mean(DF$Shannon[(DF$site=='colon')&(DF$category=='Control')])
mean(DF$InvSimpson[(DF$site=='colon')&(DF$category=='EGS')])/mean(DF$InvSimpson[(DF$site=='colon')&(DF$category=='Control')])
mean(DF$Fisher[(DF$site=='caecum')&(DF$category=='EGS')])/mean(DF$Fisher[(DF$site=='caecum')&(DF$category=='Control')])
mean(DF$Fisher[(DF$site=='colon')&(DF$category=='EGS')])/mean(DF$Fisher[(DF$site=='colon')&(DF$category=='Control')])


#now start beta diversity: Plot the PCoA using the unweighted UniFrac as distance: (these functions are not used eventually)
#Principal coordinates analysis (PCoA) derived from unweighted and weighted UniFrac and Bray-Curtis distances among samples
select_sample<-function(physeqobj,category_string){
physeq<-phyloseq(
otu_table(physeqobj)[,sample_data(physeqobj)$category==category_string],
sample_data(physeqobj)[sample_data(physeqobj)$category==category_string],
tax_table(physeqobj),
phy_tree(physeqobj))
sample_data(physeq)$category_site<-factor(sample_data(physeq)$category_site,levels = c('Control_stomach','EGS_stomach','Control_ileum','EGS_ileum','Control_caecum','EGS_caecum','Control_colon','EGS_colon','Control_faeces','EGS_faeces','Co-Grazing_faeces'))
sample_data(physeq)$site<-factor(sample_data(physeq)$site,levels = c('stomach','ileum','caecum','colon','faeces'))
return(physeq)
}
#
plot_pcoa<-function(physeqobj,wunifrac_dist,category_string,method_string){
ordination = ordinate(select_sample(physeqobj,category_string), method=method_string, distance=wunifrac_dist)
plot_ordination(select_sample(physeqobj,category_string), ordination, color="site",title=paste('Beta Plot of',category_string,'samples using',method_string,sep=' ')) + theme(aspect.ratio=1)
}
#get the graphs
dist_egs <- phyloseq::distance(select_sample(physeq1,'EGS'), method="unifrac", weighted=F)
dist_ctrl <- phyloseq::distance(select_sample(physeq1,'Control'), method="unifrac", weighted=F)
pdf('beta.pdf',height=8,width=10)
plot_pcoa(physeq1,dist_egs,'EGS','DCA')
plot_pcoa(physeq1,dist_ctrl,'Control','DCA')
plot_pcoa(physeq1,dist_egs,'EGS','CCA')
plot_pcoa(physeq1,dist_ctrl,'Control','CCA')
plot_pcoa(physeq1,dist_egs,'EGS','RDA')
plot_pcoa(physeq1,dist_ctrl,'Control','RDA')
plot_pcoa(physeq1,dist_egs,'EGS','NMDS')
plot_pcoa(physeq1,dist_ctrl,'Control','NMDS')
plot_pcoa(physeq1,dist_egs,'EGS','MDS')
plot_pcoa(physeq1,dist_ctrl,'Control','MDS')
plot_pcoa(physeq1,dist_egs,'EGS','PCoA')
plot_pcoa(physeq1,dist_ctrl,'Control','PCoA')
dev.off()
