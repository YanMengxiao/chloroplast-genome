
setwd("D:/research/ipomoea/chloroplast/alignment/snp")

############  divide group  #########
# library(ape)
# tree<-read.tree('D:/research/ipomoea/chloroplast/tree/ceg11/re_assembly/plot/tree/reroot_RAxML_bipartitions.new_ceg11_out4_BS1000')
# # ggtree(raxml) + geom_text2(aes(subset=!isTip, label=node,color="red"), hjust=-.3) + geom_tiplab(size=2)
# # ggsave("raxml_preview.pdf",width=20,height=10)node_raxml=c(119,172,203)
# node=c(172,119,203,114)
# tree<-groupClade(tree,node)
# group<-data.frame(edge=1:length(attr(tree,"group")),group=as.numeric(attr(tree,"group")))
# group<-group[1:length(tree$tip.label),]
# group$tip=tree$tip.label
# group$group<-gsub("2","B1",group$group)
# group$group<-gsub("3","B2",group$group)
# group$group<-gsub("4","T",group$group)
# group$group<-gsub("5","OUT",group$group)
# group$group<-factor(group$group)
# write.csv(group,"B1_B2_T_group.csv")
group<-read.csv("B1_B2_T_group.csv")
b1_num<-length(grep("B1",group$group))
b2_num<-length(grep("B2",group$group))
t_num<-length(grep("T",group$group))

#############  find variable site ###############
vcf<-read.csv("rob_ceg11.csv")
dim(vcf)

t_col<-match(group$tip[(grep("T",group$group))],colnames(vcf))
b1_col<-match(group$tip[(grep("B1",group$group))],colnames(vcf))
b2_col<-match(group$tip[(grep("B2",group$group))],colnames(vcf))

t<-c()
b1<-c()
b2<-c()
for (i in 1:nrow(vcf)) {
  if (sum(vcf[i,b1_col])<5 | sum(vcf[i,b1_col])>b1_num-5) {t=c(t,i)}
  if (sum(vcf[i,b2_col])<5 | sum(vcf[i,b2_col])>b2_num-5) {b1=c(b1,i)}
  if (sum(vcf[i,t_col])<5 | sum(vcf[i,t_col])>t_num-5) {b2=c(b2,i)}
}

num_filter<-sort(unique(t,b1,b2))
vcf_filter<-vcf[num_filter,]
# È¡ÖÚÊý
mode = function(x) 
{as.numeric(names(table(x))[table(x)==max(table(x))])}

site=list()
for (i in 1:nrow(vcf_filter)) {
  if (((length(mode(as.numeric(vcf_filter[i,t_col])))==1) && (length(mode(as.numeric(vcf_filter[i,b1_col])))==1)) && (length(mode(as.numeric(vcf_filter[i,b2_col])))==1))
  {if((mode((as.numeric(vcf_filter[i,t_col]))))!=(mode((as.numeric(vcf_filter[i,b1_col]))))) {site$tb1<-c(site$tb1,i)}
    if((mode((as.numeric(vcf_filter[i,t_col]))))!=(mode((as.numeric(vcf_filter[i,b2_col]))))) {site$tb2<-c(site$tb2,i)}
    if((mode((as.numeric(vcf_filter[i,b1_col]))))!=(mode((as.numeric(vcf_filter[i,b2_col]))))) {site$b1b2<-c(site$b1b2,i)}
  }
}
site_align<-lapply(site,function(x) {(vcf_filter$POS[sort(x)])[(x>17)&(x<1923)]})


################ locate site to reference ##############
# source("http://bioconductor.org/biocLite.R")
# biocLite("Biostrings") 
library("Biostrings")
library(stringr)

align_ref <- readDNAStringSet("align_xushu18.fasta")

loc_snp_site<-function(x) {
  for (i in 1:length(x)) {
    dif<-str_count(as.character(subseq(align_ref, start = 1, end = x[i])),"-")
    x[i]<-x[i]-dif
  }
  return(unique(x))}
site_ref<-lapply(site_align,loc_snp_site)

# method 2
# for (i in 1:length(site_align$tb1)) {
#   dif<-str_count(as.character(subseq(align_ref, start = 1, end = site_align$tb1[i])),"-")
#   site_align$tb1[i]<-site_align$tb1[i]-dif
# }
# for (i in 1:length(site_align$tb2)) {
#   dif<-str_count(as.character(subseq(align_ref, start = 1, end = site_align$tb2[i])),"-")
#   site_align$tb2[i]<-site_align$tb2[i]-dif
# }
# for (i in 1:length(site_align$b1b2)) {
#   dif<-str_count(as.character(subseq(align_ref, start = 1, end = site_align$b1b2[i])),"-")
#   site_align$b1b2[i]<-site_align$b1b2[i]-dif
# }
# site_ref<-site_align
# site_ref<-lapply(site_ref,unique)



###############   build gene table    ##################

# genbank<-read.csv("NCBI_xushu18.csv",col.names = c("num","info"))
# gene_info<-data.frame(range=genbank$info[grep("     gene            ",genbank$info)])
# gene_info$range<-gsub("[a-z].*([a-z]| )([0-9].*[0-9]$)","\\2",gene_info$range)
# gene_info$start<-gsub("([0-9].*[0-9])\\.\\.([0-9].*[0-9])","\\1",gene_info$range)
# gene_info$end<-gsub("([0-9].*[0-9])\\.\\.([0-9].*[0-9])","\\2",gene_info$range)
# gene_info$name<-genbank$info[grep("     gene            ",genbank$info)+1]
# gene_info$name<-gsub("/gene=","",gene_info$name)
# gene_info$name<-gsub(" ","",gene_info$name)
# write.csv(gene_info,"gene_info.csv")

gene_info<-read.csv("gene_info.csv",as.is =T)
loci<-c()
loci[grep("-",gene_info$name)]<-"intergenic"
loci[!(1:nrow(gene_info))%in%grep("-",gene_info$name)]<-"CDS"
gene_info$loci<-loci

snp2gene<-function(x) {
  tb1<-list()
for (snp in x) {
  for (n in 1:nrow(gene_info)) {
    if (snp>gene_info$start[n] & snp<gene_info$end[n]) 
    {tb1$name<-append(tb1$name,gene_info$name[n])
    tb1$site<-append(tb1$site,snp)
    tb1$region<-append(tb1$region,gene_info$region[n])
    tb1$loci<-append(tb1$loci,gene_info$loci[n])}
    else if (snp>gene_info$end[n] & snp<gene_info$start[n+1]) 
    {tb1$name<-append(tb1$name,paste(gene_info$name[n],"-",gene_info$name[n+1]))
    tb1$site<-append(tb1$site,snp)
    tb1$region<-append(tb1$region,gene_info$region[n])
    tb1$loci<-append(tb1$loci,gene_info$loci[n])}
  }
}
tb1<-data.frame(site=tb1$site,name=tb1$name,region=tb1$region,loci=tb1$loci)
return(tb1)
}

snp_gene<-lapply(site_ref,snp2gene)
list2matrix<-function(x) {
  table<-rbind(rbind(x[[1]],x[[2]]),x[[3]]) 
  table$group<-c(rep(names(x[1]),nrow(x[[1]])),rep(names(x[2]),nrow(x[[2]])),rep(names(x[3]),nrow(x[[3]]))) 
  return(table)}
# sum_gene<-list2matrix(snp_gene)
# tb1<-snp2gene(site_ref$tb1)
# tb2<-snp2gene(site_ref$tb2)
# b1b2<-snp2gene(site_ref$tb3)
write.csv(snp_gene$tb1,"tb1_gene.csv")
write.csv(snp_gene$tb2,"tb2_gene.csv")
write.csv(snp_gene$b1b2,"b1b2_gene.csv")


sum_gene_name<-function(x) {
sum_tb1<-data.frame(sort(table(x$name),decreasing =T))
for (i in 1:nrow(sum_tb1)) {
  sum_tb1$site[i]<-paste(x$site[grep(sum_tb1$Var1[i],x$name)],collapse = ";")
}
sum_tb1$region<-x$region[match(sum_tb1$Var1,x$name)]
sum_tb1$loci<-x$loci[match(sum_tb1$Var1,x$name)]
return(sum_tb1)
}
summary_gene<-lapply(snp_gene, sum_gene_name)
sum_gene<-rbind(rbind(summary_gene$tb1,summary_gene$tb2),summary_gene$b1b2)
sum_gene$group<-c(rep("TB1",nrow(summary_gene$tb1)),rep("TB2",nrow(summary_gene$tb2)),rep("B1B2",nrow(summary_gene$b1b2)))
sum_gene$CS<-sum_gene$region
sum_gene$CS<-gsub("L|A|B","",sum_gene$CS)
sum_gene$CS<-gsub("SSC","SC",sum_gene$CS)
write.csv(sum_gene,"site_gene_summary2.csv")

sum_gene<-read.csv("site_gene_summary2.csv")


summary_region<-lapply(summary_gene,function(x){aggregate(x$Freq, by= list(x$region), FUN = sum)})
sum_region<-list2matrix(summary_region)
colnames(sum_region)<-c("region","Freq","group")
summary_loci<-lapply(summary_gene,function(x){aggregate(x$Freq, by= list(x$loci), FUN = sum)})
sum_loci<-list2matrix(summary_loci)
colnames(sum_loci)<-c("region","Freq","group")

write.csv(sum_region,"summary_region.csv")
sum_region<-read.csv("summary_region.csv")

#############   plot  ################

library(ggplot2)
summary<-read.csv("site_gene_summary.csv")
ggplot(data=summary,aes(x=group,y=gene)) + geom_point(size=summary$Freq,color="blue")

ggplot(data=sum_region) +
  geom_bar(aes(x=group,y=Freq,fill=region),stat="identity",position="dodge",width=0.5)+
  labs(y="number",fill="region")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_line(colour = NA),panel.grid.minor = element_line(colour = NA))
ggsave("sum_region.pdf")

ggplot(data=sum_loci) +
  geom_bar(aes(x=group,y=Freq,fill=region),stat="identity",position="dodge",width=0.5)+
  labs(y="number",fill="region")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_line(colour = NA),panel.grid.minor = element_line(colour = NA))
ggsave("sum_loci.pdf")

ggplot(data=sum_gene) +
  geom_bar(aes(x=group,y=Freq,fill=CS),stat="identity",position="dodge",width=0.5)+
  labs(y="number",fill="region")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_line(colour = NA),panel.grid.minor = element_line(colour = NA))
ggsave("sum_single_copy.pdf")

table <- read.csv(file ='resolution.csv')
table$Marker=factor(table$Marker, levels=c("A","M","T","Y","AM","AT","AY","MT","MY","TY","AMT","AMY","ATY","MTY","AMTY"))
table$method=factor(table$method, levels=c("BLOG","Tree","Distance","BLAST"))
discrimination_rate_plot <-ggplot(table) + 
  geom_bar(aes(x=Marker,y=value,fill=method),stat="identity",position="dodge")+ 
  scale_fill_brewer(palette="gray") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  labs(y="Discrimination rate (%)")

# sum_gene_region<-function(x) {
#   sum_tb1<-data.frame(sort(table(x$region),decreasing =T))
#   return(sum_tb1)
# }
# summary_gene_region<-lapply(snp_gene, sum_gene_region)
# sum_gene_region<-rbind(rbind(summary_gene_region$tb1,summary_gene_region$tb2),summary_gene_region$b1b2)
# sum_gene_region$group<-c(rep("TB1",4),rep("TB2",4),rep("B1B2",4))
# 
# sum_gene_loci<-function(x) {
#   sum_tb1<-data.frame(sort(table(x$loci),decreasing =T))
#   return(sum_tb1)
# }
# summary_gene_loci<-lapply(snp_gene, sum_gene_loci)
# sum_gene_loci<-rbind(rbind(summary_gene_loci$tb1,summary_gene_loci$tb2),summary_gene_loci$b1b2)
# sum_gene_loci$group<-c(rep("TB1",2),rep("TB2",2),rep("B1B2",2))


# FAIL!!!!!!!!!!!!
# sum_gene<-function(x,y) {
#   sum_tb1<-data.frame(sort(table(x$y),decreasing =T))
#   for (i in 1:nrow(sum_tb1)) {
#     sum_tb1$site[i]<-paste(x$site[grep(as.character(sum_tb1$Var1[i]),x$name)],collapse = ";")
#   }
#   sum_tb1$region<-x$region[match(sum_tb1$Var1,x$name)]
#   return(sum_tb1)
# }
# summary_gene<-lapply(snp_gene, sum_gene, y=c("name","region","loci"))



# sum_tb1<-data.frame(sort(table(tb1$name),decreasing =T))
# for (i in 1:nrow(sum_tb1)) {
#   sum_tb1$site[i]<-paste(tb1$site[grep(sum_tb1$Var1[i],tb1$name)],collapse = ";")
# }
# sum_tb1$region<-tb1$region[match(sum_tb1$Var1,tb1$name)]
# 
# sum_tb2<-data.frame(sort(table(tb2$name),decreasing =T))
# for (i in 1:nrow(sum_tb2)) {
#   sum_tb2$site[i]<-paste(tb2$site[grep(sum_tb2$Var1[i],tb2$name)],collapse = ";")
# }
# sum_tb2$region<-tb2$region[match(sum_tb2$Var1,tb2$name)]
# 
# sum_b1b2<-data.frame(sort(table(b1b2$name),decreasing =T))
# for (i in 1:nrow(sum_b1b2)) {
#   sum_b1b2$site[i]<-paste(b1b2$site[grep(sum_b1b2$Var1[i],b1b2$name)],collapse = ";")
# }
# sum_b1b2$region<-b1b2$region[match(sum_b1b2$Var1,b1b2$name)]

# sum_reg<-data.frame(sort(table(tb1$name),decreasing =T))
# 
# summary<-rbind(rbind(sum_tb1,sum_tb2),sum_b1b2)
# summary$group<-c(rep("TB1",nrow(sum_tb1)),rep("TB2",nrow(sum_tb2)),rep("B1B2",nrow(sum_b1b2)))
# colnames(summary)[1]<-"gene"
# # summary<-merge(merge(sum_tb1,sum_tb2,by="Var1",all=T),sum_b1b2,by="Var1",all=T)
# # colnames(summary)<-c("gene","Freq_TB1","Site_TB1","Freq_TB2","Site_TB2","Freq_B1B2","Site_B1B2")
# write.csv(summary,"site_gene_summary2.csv")
# summary<-read.csv("site_gene_summary.csv")


############### reorder seq by groups  #################

seq<- data.frame(name = names(fastaFile), sequence = paste(fastaFile))
seq$group<-group$group[match(seq$name,group$tip)]
seq<-seq[order(seq$group),]


sort.seq<-paste(">",seq$name,sep="")
sort.seq<-paste(sort.seq,seq$sequence,sep="###")
write.csv(sort.seq,"sort.seq.csv") 
