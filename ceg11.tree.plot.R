# ipomoea
## gblock10
library(ggtree)
library(ape)
library(phytools)
library(cowplot)
library(phangorn)

# ipomoea gblock5
setwd('D:/research/ipomoea/chloroplast/tree/ceg11/re_assembly/plot')


#######  raxml_GTR+I+G_BS1000 ####### 
raxml<-read.nexus("tree/root_raxml_bipartitions.new_ceg11_out4_bs1000.tre")
raxml<-drop.tip(raxml,c("I_leucantha_CIP460002","I_triloba_CIP460554"))
raxml$node.label<-gsub("\'","",raxml$node.label)
# raxml<-root(raxml,"I_cynanchifolia_JRIW_27533")
# raxml<-midpoint(raxml)
ggtree(raxml) + geom_nodelab()
## manual check group node number
ggtree(raxml) + geom_text2(aes(subset=!isTip, label=node,color="red"), hjust=-.3) + geom_tiplab(size=2)
ggsave("raxml_preview.pdf",width=20,height=10)
## viewClade(ggtree(tree)+geom_tiplab(size=2), node=250)
## label new tips
node_raxml=c(115,189,168)
raxml<- groupClade(raxml, node_raxml)
tip_col_raxml=data.frame(edge=1:219,group2=as.numeric(attr(raxml,"group")),group=factor(as.numeric(attr(raxml,"group"))))
tip_col_raxml$group2[grep("4X",raxml$tip.label)]<-5
tip_col_raxml$group2[grep("Taizhong6|CIP420065|PI561544|CIP460814|KP212149",raxml$tip.label)]<-6
tip_col_raxml$group=factor(tip_col_raxml$group2-1)
plot_raxml<-ggtree(raxml,aes(color=attr(raxml,"group")))+
  scale_color_manual(values=c("gray50","#F8766D", "#00B6EB", "#FAA200","#00C094","#A58AFF"))
tip_raxml<-plot_raxml %<+% tip_col_raxml + geom_tiplab(aes(color=factor(tip_col_raxml$group)),size=2,align=T,linesize=0.5)
## select bootstrap > 70
dr <- tip_raxml$data
dr <- dr[!dr$isTip,]
dr$label <- as.numeric(dr$label)
dr <- dr[dr$label > 70,]
final_raxml<- tip_raxml + geom_text(data=dr, aes(label=label),size = 3, col= "black",hjust=-.3)+ 
  xlim(0, 0.0015)+ labs(title="raxml")
ggsave("raxml_ceg11_new_gblock5.pdf",final_raxml,width=20,height=10)

####### iqtree_GTR+I+G_BS1000  #######
iqtree<-read.nexus("tree/root_gblock5_bt_out4_ceg11_new_assembly.fasta.treefile.tre")
iqtree<-drop.tip(iqtree,c("I_leucantha_CIP460002","I_triloba_CIP460554"))
iqtree$node.label<-gsub("\'","",iqtree$node.label)
# manual check group node number
ggtree(iqtree) + geom_text2(aes(subset=!isTip, label=node,color="red"), hjust=-.3) + geom_tiplab(size=2)
ggsave("iqtree_preview.pdf",width=20,height=10)
# viewClade(ggtree(tree)+geom_tiplab(size=2), node=250)
# group tree and color subtree accoring to groups
## label new tips
node_iqtree=c(116,189,168)
iqtree<- groupClade(iqtree, node_iqtree)
tip_col_iqtree=data.frame(edge=1:219,group2=as.numeric(attr(iqtree,"group")),group=factor(as.numeric(attr(iqtree,"group"))))
tip_col_iqtree$group2[grep("4X",iqtree$tip.label)]<-5
tip_col_iqtree$group2[grep("Taizhong6|CIP420065|PI561544|CIP460814|KP212149",iqtree$tip.label)]<-6
tip_col_iqtree$group=factor(tip_col_iqtree$group2-1)
plot_iqtree<-ggtree(iqtree,aes(color=attr(iqtree,"group")))+
  scale_color_manual(values=c("gray50","#F8766D", "#00B6EB", "#FAA200","#00C094","#A58AFF"))
tip_iqtree<-plot_iqtree %<+% tip_col_iqtree + geom_tiplab(aes(color=factor(tip_col_iqtree$group)),size=2,align=T,linesize=0.5)
# select bootstrap > 70
di <- tip_iqtree$data
di <- di[!di$isTip,]
di$label <- as.numeric(di$label)
di <- di[di$label > 70,]
final_iqtree<- tip_iqtree + geom_text(data=di, aes(label=label),size = 3, col= "black",hjust=-.3)+ 
  xlim(0, 0.0015)+ labs(title="iqtree")
ggsave("iqtree_ceg11_gblock5.pdf",final_iqtree,width=20,height=10)

#######  NJ_BS1000  ####### 
nj<-read.nexus("tree/root_nj_1000_gblock5_bt_out4_ceg11_new_assembly-18936.nwk")
nj<-drop.tip(nj,c("I_cordatotriloba_JRIW_27733","I_cynanchifolia_JRIW_27533"))
nj$node.label<-gsub("\'","",nj$node.label)
# manual check group node number
ggtree(nj) + geom_text2(aes(subset=!isTip, label=node,color="red"), hjust=-.3) + geom_tiplab(size=2)
ggsave("nj_preview.pdf",width=20,height=10)
# viewClade(ggtree(tree)+geom_tiplab(size=2), node=250)
## label new tips
node_nj=c(168,137,115)
nj<- groupClade(nj, node_nj)
tip_col_nj=data.frame(edge=1:219,group2=as.numeric(attr(nj,"group")),group=factor(as.numeric(attr(nj,"group"))))
tip_col_nj$group2[grep("4X",nj$tip.label)]<-5
tip_col_nj$group2[grep("Taizhong6|CIP420065|PI561544|CIP460814|KP212149",nj$tip.label)]<-6
tip_col_nj$group=factor(tip_col_nj$group2-1)
plot_nj<-ggtree(nj,aes(color=attr(nj,"group")))+
  scale_color_manual(values=c("gray50","#F8766D", "#00B6EB", "#FAA200","#00C094","#A58AFF"))
tip_nj<-plot_nj %<+% tip_col_nj + geom_tiplab(aes(color=factor(tip_col_nj$group)),size=2,align=T,linesize=0.5)
# select bootstrap > 70
dnj <- tip_nj$data
dnj <- dnj[!dnj$isTip,]
dnj$label <- as.numeric(dnj$label)
dnj <- dnj[dnj$label > 0.7,]
final_nj<- tip_nj + geom_text(data=dnj, aes(label=round(label,2)*100),size = 3, col= "black",hjust=-.3)+ 
  xlim(0, 0.0015)+ labs(title="nj")
ggsave("nj_robert_ipomoea_gblock5.pdf",final_nj,width=20,height=10)


#######  MP_BS1000  ####### 
mp<-read.nexus("tree/root_mp_1000_gblock5_BT_out4_ceg11_new_assembly-26665.nwk")
mp<-drop.tip(mp,c("I_leucantha_CIP460002","I_triloba_CIP460554"))
mp$node.label<-gsub("\'","",mp$node.label)
# manual check group node number
ggtree(mp) + geom_text2(aes(subset=!isTip, label=node,color="red"), hjust=-.3) + geom_tiplab(size=2)
ggsave("mp_preview.pdf",width=20,height=10)
# viewClade(ggtree(tree)+geom_tiplab(size=2), node=250)
## label new tips
node_mp=c(137,189,115)
mp<- groupClade(mp, node_mp)
tip_col_mp=data.frame(edge=1:219,group2=as.numeric(attr(mp,"group")),group3=attr(mp,"group"))
tip_col_mp$group2[grep("4X",mp$tip.label)]<-5
tip_col_mp$group2[grep("Taizhong6|CIP420065|PI561544|CIP460814|KP212149",mp$tip.label)]<-6
tip_col_mp$group=factor(tip_col_mp$group2-1)
plot_mp<-ggtree(mp,aes(color=attr(mp,"group")))+
  scale_color_manual(values=c("gray50","#F8766D", "#00B6EB", "#FAA200","#00C094","#A58AFF"))
tip_mp<-plot_mp %<+% tip_col_mp + geom_tiplab(aes(color=factor(tip_col_mp$group)),size=2,align=T,linesize=0.5)
# select bootstrap > 70
dmp <- tip_mp$data
dmp <- dmp[!dmp$isTip,]
dmp$label <- as.numeric(dmp$label)
dmp <- dmp[dmp$label > 0.7,]
final_mp<- tip_mp + geom_text(data=dmp, aes(label=round(label,2)*100),size = 3, col= "black",hjust=-.3)+
  xlim(0, 300)+ labs(title="mp")
ggsave("mp_robert_BT_gblock5.pdf",final_mp,width=20,height=10)

ml<-plot_grid(final_raxml,final_iqtree,ncol=2)
ggsave("ML_ceg11_gblock5.pdf",ml,width=20,height=10)
nj_mp<-plot_grid(final_nj,final_mp,ncol=2)
ggsave("NJ_MP_ceg11_gblock5.pdf",nj_mp,width=20,height=10)
four<-plot_grid(final_raxml,final_iqtree,final_nj,final_mp,ncol=4)
ggsave("four_ceg11_new_gblock5.pdf",four,width=20,height=10)

