library(readr)
library(plyr)
library(forcats)
library(ggplot2)
library(ggsci)
library(dplyr)
library(ggalluvial)
library(RColorBrewer)
library(reshape2)
library(vegan)
library(tidyr)
library(ggsignif)
library(scales)
library(ape)
library(tidyverse)
library(patchwork)
library(ggpmisc)
library(ggpubr)


# 1. Abundance ------------------------------------------------------------
filtered <- read.csv("/mnt/raid8/sunchuqing/2022_OtherAna/Filtered.CHGV/filtered.csv",quote="")
kept<-subset(filtered,Annotation=="≥2 software" |(Annotation=="≥1 software" & checkv_quality%in%c("High-quality","Complete")))$Contig %>% unique 
positive <- read.csv("/mnt/raid8/sunchuqing/2022_OtherAna/Filtered.CHGV/positive.csv",quote="",header=F)
circular <- c(subset(positive,V1=="cir")[,2],subset(filtered,completeness==100)[,1])
source("00_Scripts/Abundance.R")
path<-"/mnt/raid8/sunchuqing/2020_Human_phage/08_Abundance_calculation"
sa<-read.table("/mnt/raid8/sunchuqing/2020_Human_phage/Dir_file.list/file.list.0818.keep",quote="",header = F)$V1

ge<-plot_VCabundance(sa,10,VC[,c(2,1)],kept)
All_abunt<-ge[[1]]
Abunt_top10<-ge[[2]]
sp<-plot_SPabundance(sa,10,kept)
All_abunt_sp<-sp[[1]]
Abunt_top10_sp<-sp[[2]]


# 2. Prevalence -----------------------------------------------------------

crass.list<-read.table("/mnt/raid8/sunchuqing/2021_OtherAna/All_contig_anno/Terminase/HGV.crass.guba/crass2.list",sep=",",header = F)
guba.list<-read.table("/mnt/raid8/sunchuqing/2021_OtherAna/All_contig_anno/Terminase/HGV.crass.guba/guba.list",sep=",",header = F)
lifestyle<-read_csv("/mnt/raid8/sunchuqing/2021_OtherAna/VC_top10_annotation/All.virulent.temperate.0729.csv",col_names = F)
lifestyle$X1<-gsub("_length.*$","",lifestyle$X1)
lifestyle <- lifestyle  %>% group_by(X1) %>% filter(X3 == max(X3)) %>% filter(X2 == max(X2))

Over.prevalence.cutoff<-data.frame(
  "cutoff"=NA,
  "VP-crAss"=NA,
  "VP-crAss-vir"=NA,
  "VC-crAss"=NA,
  "VP-guba"=NA,
  "VP-guba-vir"=NA,
  "VC-guba"=NA
)

for(cutoff in c(1e-8,1e-7,1e-6,1e-5,1e-4)){
  prevalence<-data.frame(table(subset(All_abunt_sp,Abundance>=cutoff)[,"id"])/length(unique(All_abunt_sp$`sample id`)))
  colnames(prevalence)<-c("id","Prevalence")
  pp.te<-prevalence
  pp.te$db<-"no"
  cr<-max(subset(pp.te,id %in%crass.list$V1)$Prevalence)
  gu<-max(subset(pp.te,id %in%guba.list$V1)$Prevalence)
  crass2<-subset(pp.te,Prevalence>cr & !id%in% c(crass.list$V1,guba.list$V1) )$id %>% unique %>% length
  guba2<-subset(pp.te,Prevalence>gu & !id%in% c(crass.list$V1,guba.list$V1) )$id %>% unique %>% length
  crass<-subset(pp.te,Prevalence>cr & id %in% subset(lifestyle,X4 %in% c("virulent","uncertain_virulent"))$X1 & (!id%in% c(crass.list$V1,guba.list$V1)))$id %>% unique %>% length
  guba<-subset(pp.te,Prevalence>gu & id %in% subset(lifestyle,X4 %in% c("virulent","uncertain_virulent"))$X1 &!id%in% c(crass.list$V1,guba.list$V1) )$id %>% unique %>% length
  prevalence.vc<-data.frame(table(subset(All_abunt,Abundance>=cutoff)[,"id"])/length(unique(All_abunt$`sample id`)))
  colnames(prevalence.vc)<-c("id","Prevalence")
  pp.vc<-prevalence.vc
  pp.vc$"Group"<-"-"
  cr.vc<-subset(VC,Contig %in%crass.list$V1)$VC %>% unique
  gu.vc<-subset(VC,Contig %in%guba.list$V1)$VC%>% unique
  
  pp.vc$Group[which(pp.vc$id %in% cr.vc)]<-'crAss'
  pp.vc$Group[which(pp.vc$id %in% gu.vc)]<-'guba'
  cr.vc<-max(subset(pp.vc,Group=="crAss")$Prevalence)
  gu.vc<-max(subset(pp.vc,Group=="guba")$Prevalence)
  
  crass.vc<-subset(pp.vc,Prevalence>cr.vc & ! Group %in% c("crAss","Guba"))$id %>% unique %>% length
  guba.vc2<-subset(pp.vc,Prevalence>gu.vc & ! Group %in% c("crAss","Guba"))$id %>% unique %>% length
  
  temp<-c(cutoff,crass2,crass,crass.vc,guba2,guba,guba.vc2)
  Over.prevalence.cutoff<-rbind(Over.prevalence.cutoff,temp)
}
na.omit(Over.prevalence.cutoff)
Over.prevalence.cutoff<-na.omit(Over.prevalence.cutoff)
Over.prevalence.cutoff



# 3. VC -------------------------------------------------------------------

mcl.input <- read.csv(file="/mnt/raid8/sunchuqing/2022_OtherAna/Filtered.CHGV/01_VC/mcl.input",
                      sep="\t",quote="",
                      header = F)
mcl.input$V1 <- gsub("_length.*$","",mcl.input$V1)
mcl.input$V2 <- gsub("_length.*$","",mcl.input$V2)
mcl.input <- subset(mcl.input,V1 %in% kept & V2 %in% kept)
write.table(mcl.input,
            file="/mnt/raid8/sunchuqing/2022_OtherAna/Filtered.CHGV/01_VC/mcl.input.2",
            sep="\t",quote=F,row.names = F,col.names = F)
write.table(kept,
            file="/mnt/raid8/sunchuqing/2022_OtherAna/Filtered.CHGV/01_VC/all.contig.list",
            sep="\t",quote=F,row.names = F,col.names = F)
VC <- read.csv(file="/mnt/raid8/sunchuqing/2022_OtherAna/Filtered.CHGV/01_VC/VC_sin.csv",
               sep=",",quote="",
               header = F)
colnames(VC)<-c("VC","Contig")

VC$Type <- "Normal"
VC$Type[which(VC$Contig %in% crass.list$V1)] <- "crAssphage"
VC$Type[which(VC$Contig %in% guba.list$V1)] <- "Gubaphage"
VC_C <- table((VC[,c("Contig","VC")] %>% unique)$VC) %>% data.frame()
colnames(VC_C)<-c("VC","Diversity")
cols <- rev(brewer.pal(9, 'PuBu'))

VC_C %>% 
  mutate(VC=fct_reorder(VC,desc(Diversity))) %>% 
  top_n(50,Diversity) %>%
  ggplot(aes(x=VC,y=Diversity,fill=Diversity))+
  geom_line(aes(group=1))+
  geom_point(aes(size=Diversity),shape=21)+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  # labs(x="Nr. of host")+
  labs(size="Nr. of host",fill="")+
  scale_fill_gradientn(colours = rev(cols))+
  scale_y_continuous(breaks = seq(0,32,2))
# 4. Viral Clade ----------------------------------------------------------

VCL.genename<-read.csv("/mnt/raid8/sunchuqing/2021_OtherAna/All_contig_anno/Terminase/CL.70.csv",header = F)
genename<-read.table("/mnt/raid8/sunchuqing/2021_Virome_base_modification/MTase.REase.ana/genename.s.re.csv",sep=",",header = F)
genename<-subset(genename,V2!="locus_tag")
genename$V1<-gsub("_length.*","",genename$V1)
genename<-genename %>% unique

VCL<-merge(VCL.genename,genename,by="V2")
VCL<-VCL[,c(2,3)] %>% unique
colnames(VCL)<-c("CL","Contig")
VCL <- subset(VCL,Contig %in% VC$Contig)
VCL$Type <- "Normal"
VCL$Type[which(VCL$Contig %in% crass.list$V1)] <- "crAssphage"
VCL$Type[which(VCL$Contig %in% guba.list$V1)] <- "Gubaphage"
VCL <- VCL %>% unique()
# VCL[,c(1,3)] %>% unique() %>% View
# VCL.rank.1 <- VCL.rank
VCL.rank <- aggregate(VCL$Contig,by=list(VCL$CL),FUN=length)
VCL.rank$New.CL[order(VCL.rank$x,decreasing = T)] <- paste("CL_",1:691,sep = "")
VCL <- merge(VCL,VCL.rank,by.x="CL",by.y="Group.1")
VCL <- VCL[,c(2,3,4,5)]
colnames(VCL)[c(3,4)] <-c("Size","CL") 
VCL.VC<- merge(VC[,c("Contig","VC")],VCL,by="Contig")
subset(VC,VC %in% VCL.VC$VC)$Contig %>% unique %>% length()
VCL.VC <- VCL.VC[,c("CL","VC","Size")] %>% unique
VCL.VC <- merge(VCL.VC,VC[,c("Contig","VC")],by="VC")
# 6. cf. NGS --------------------------------------------------------------
VC <- merge(VC,lifestyle[,c(1,2,4)],by.x="Contig",by.y="X1") %>% unique()
colnames(VC)[c(4,5)] <- c("Length","Lifestyle") 
VC$Assembly <- "NGS"
VC$Assembly[grep("NODE",VC$Contig)] <- "Hybrid"
VC$Assembly[grep("opera",VC$Contig)] <- "Hybrid"
VC$Assembly[grep("pilon",VC$Contig)] <- "TGS"
VC$Length <- as.numeric(VC$Length)
VC$Length[which(VC$Contig=="DC05_NODE_305")] <- 6724
VC$Length[which(VC$Contig=="SK43_NODE_80")] <- 43135
VC %>% 
  mutate(Assembly=factor(Assembly,levels = c("NGS","Hybrid","TGS"))) %>% 
  mutate(Complete=factor(Complete,levels = c("Linear","Complete","Circular"))) %>% 
  
  ggplot(aes(x=Assembly,y=Length,color=Complete))+
  geom_boxplot(width=0.5,outlier.size = 0.6)+
  scale_y_log10()+
  theme_bw()+
  labs(y="Genome Length",x="Assembly Strategy",color="")+
  scale_color_npg()
#


# 6.5 Circular/Complete  --------------------------------------------------
VC.bedtools <- VC[,c(1,4)] %>% unique()
colnames(VC.bedtools) <- c("Contig","End")
VC.bedtools$Start <- VC.bedtools$End-1
VC.bedtools.end <- VC.bedtools[,c(1,3,2)]
VC.bedtools.start <- data.frame(Contig=VC.bedtools$Contig,
                                Start=1,
                                End=2)
write.table(VC.bedtools.start,"/mnt/raid8/sunchuqing/2022_OtherAna/Filtered.CHGV/04_mapping/bedtools.start",col.names=F,sep="\t",quote = F,row.names = F)
write.table(VC.bedtools.end,"/mnt/raid8/sunchuqing/2022_OtherAna/Filtered.CHGV/04_mapping/bedtools.end",col.names=F,sep="\t",quote = F,row.names = F)



VC$Complete <- "Others"
circular2 <- read.table("/mnt/raid8/sunchuqing/2022_OtherAna/Filtered.CHGV/Mapping.circular",sep = ",",quote="") %>% unique()
circular2 <- (table(circular2$V2) %>% as.data.frame() %>% subset(Freq>1))[,1]
subset(VC,Contig %in% circular2 | Complete  %in% c("Complete","Circular"))$Contig %>% length()
subset(VC,Contig %in% circular2 &(! Complete  %in% c("Complete","Circular")))$Contig %>% length()

circular2 %>% length

VC$Complete[which(VC$Contig %in% subset(filtered,completeness==100)[,1])] <- "Complete"
VC$Complete[which(VC$Contig %in% subset(positive,V1=="cir")[,2])] <- "Circular"
VC$Complete[which(VC$Contig %in% circular2)] <- "Circular"

VC$Annotation <- "≥1 software"
VC$Annotation[which(VC$Contig %in% subset(filtered,Annotation=="≥2 software")[,1])]  <- "≥2 software"
links<-aggregate(unique(VC[,c("Contig","Complete","Assembly","Annotation")])$Complete,
                 list(
                   unique(VC[,c("Contig","Complete","Assembly","Annotation")])$Complete,
                   # unique(VC[,c("Contig","Complete","Assembly","Annotation")])$Assembly,
                   unique(VC[,c("Contig","Complete","Assembly","Annotation")])$Annotation),FUN=length) 
# names(links)<-c("Assembly","Complete","Annotation","value")
names(links)<-c("Complete","Annotation","value")

links %>% 
  # mutate(Assembly=factor(Assembly,levels = c("NGS","Hybrid","TGS"))) %>% 
  mutate(Complete=factor(Complete,levels = c("Others","Circular","Complete"))) %>% 
  ggplot(aes( axis1 = Annotation, axis2=Complete, y= value,)) +
  scale_x_discrete(limits = c("Annotaion","Complete"), expand = c(.1, .05)) +
  # ggplot(aes(axis1 = Assembly, axis2 = Complete, axis3=Annotation, y= value,)) +
  # scale_x_discrete(limits = c("Assembly", "Complete","Annotaion"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Complete)) +
  geom_stratum() + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  theme_bw() +
  scale_fill_igv()+
  theme(plot.margin=unit(rep(1,4),'lines'))+
  theme(strip.background.x = element_rect(fill = "white")) 
VC$Complete <- "All"
VC$Complete[which(VC$Contig %in% c(subset(filtered,completeness==100)[,1],subset(positive,V1=="cir")[,2]))] <- "Complete"
tmp<- VC[which(VC$Contig %in% c(subset(filtered,completeness==100)[,1],subset(positive,V1=="cir")[,2])),] 
VC$Complete <- "All"
tmp2 <- VC
tmp2$Assembly <- "CHGV(Combined)"
tmp3 <- tmp
tmp3$Assembly <- "CHGV(Combined)"
# rbind(VC,tmp) %>% 
#   mutate(Assembly=factor(Assembly,levels = c("NGS","Hybrid","TGS"))) %>% 
#   ggplot(aes(x="Assembly",y=..count..,fill=Complete))+
#   geom_bar(width=0.7,position = "fill",color="white")+
#   theme_void()+
#   # theme(axis.x.test = element_blank()) +
#   labs(y="",fill="")+
#   facet_grid(~Assembly)+
#   coord_polar(theta = 'y')+
#   scale_y_continuous(labels = scales::percent)+
#   scale_fill_simpsons()
le <- rbind(subset(VC,Assembly %in% c("NGS","TGS")),subset(tmp,Assembly %in% c("NGS","TGS")),tmp2,tmp3)[,c("Assembly","Complete","Length")] 
GVD.meta$Complete <- "All"
gvd.tmp <- subset(GVD.meta,completeness==100)
gvd.tmp$Complete <- "Complete"
gvd <- rbind(GVD.meta,gvd.tmp)[,c("Length","Complete")]
colnames(gvd)[1] <- "Length"
gvd$"Assembly" <- "GVD"
GPD.checkv <- read.csv("/mnt/raid8/sunchuqing/Database/GPD/GPD_metadata.tsv",sep="\t")
GPD.checkv$Complete <- "All"
gpd.tmp <- subset(GPD.checkv,checkV_completion==100)
gpd.tmp$Complete <- "Complete"
gpd <- rbind(GPD.checkv,gpd.tmp)[,c("Size","Complete")]
colnames(gpd)[1] <- "Length"
gpd$"Assembly" <- "GPD"
rbind(le,gvd,gpd)%>% 
  mutate(Assembly=factor(Assembly,levels = c("NGS","TGS","CHGV(Combined)","GVD","GPD"))) %>% 
  mutate(Complete=factor(Complete,levels = c("All","Complete"))) %>% 
  
  ggplot(aes(x=Assembly,y=Length,color=Complete))+
  geom_boxplot(width=0.5,outlier.size = 0.6)+
  scale_y_log10()+
  theme_bw()+
  labs(y="Genome Length",x="Assembly Strategy",color="")+
  scale_color_npg()

rbind(le,gvd,gpd)%>% 
  subset(Length>=10000) %>% 
  mutate(Assembly=factor(Assembly,levels = c("NGS","TGS","CHGV(Combined)","GVD","GPD"))) %>% 
  mutate(Complete=factor(Complete,levels = c("All","Complete"))) %>% 
  
  ggplot(aes(x=Assembly,y=Length,color=Complete))+
  geom_signif(comparisons = list(c("CHGV(Combined)","GPD")),
              color="Black",
              step_increase = 0.1,
              map_signif_level = T,test = t.test)+
  geom_boxplot(width=0.5,outlier.size = 0.6)+
  scale_y_log10()+
  theme_bw()+
  labs(y="Genome Length(>10k)",x="Assembly Strategy",color="")+
  scale_color_npg()
VC<- VC %>% 
  mutate(Assembly=factor(Assembly,levels = c("NGS","Hybrid","TGS"))) 
per <- aggregate(VC$Complete,by=list(VC$Assembly,VC$Complete),FUN=length)
colnames(per) <- c("Assembly","Complete","num")
per <- rbind(per,data.frame(Assembly="GPD",Complete=c("Complete","Linear"),num=c(17743,142810-17743)))
per%>% 
  ggplot(aes(x="Assembly",y=num,fill=Complete))+
  geom_col(width=0.7,position = "fill",color="white")+
  theme_void()+
  # theme(axis.x.test = element_blank()) +
  labs(y="",fill="")+
  facet_wrap(~Assembly,ncol=2)+
  coord_polar(theta = 'y')+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_simpsons()

# 7. Shannon --------------------------------------------------------------
vcutoff<-1e-6
bcutoff<-1e-6
abb<-subset(All_abunt,Abundance>=vcutoff)[,c(1,2,3)]
VC_abundance_cast<-spread(abb,`sample id`,Abundance)
VC_abundance_cast[is.na(VC_abundance_cast)]<-0
rownames(VC_abundance_cast)<-VC_abundance_cast[,1]
VC_abundance_cast<-VC_abundance_cast[,-1]
# VC_abundance_cast<-VC_abundance_cast[which(rowSums(VC_abundance_cast) > (vcutoff*10)),]
# VC_abundance_cast<-t(VC_abundance_cast)
# VC_abundance_cast<-VC_abundance_cast[which(rowSums(VC_abundance_cast) > (vcutoff*10)),]
# VC_abundance_cast<-t(VC_abundance_cast)
VC_abundance_cast<-sweep(VC_abundance_cast,2,colSums(VC_abundance_cast),`/`)
VC_abundance_cast<-t(VC_abundance_cast) %>% melt()
names(VC_abundance_cast)<-c("Sample","id","Abundance")
abb<-VC_abundance_cast
prevalence<-data.frame(table(subset(abb,Abundance>=vcutoff)[,"id"])/length(unique(abb$Sample)))

colnames(prevalence)<-c("id","Prevalence")
bacprevalence2<-read_csv("/mnt/raid8/sunchuqing/2020_Human_phage/0001_BacMetagenomic/06_Abundance/genus_abundance_profile.csv")
bacprevalence2<-data.frame(bacprevalence2)
bp2<-bacprevalence2[,-1]
rownames(bp2)<-bacprevalence2$taxonomy

# bp2<-bp2[which(rowSums(bp2) > (bcutoff*10)),]
# bp2<-t(bp2)
# bp2<-bp2[which(rowSums(bp2) > (bcutoff*10)),]
# bp2<-t(bp2)
bp2<-bp2 %>%t() %>%melt() 
names(bp2)<-c("Sample","id","Abundance")
bp2$Abundance <- bp2$Abundance/100

bp2<-data.frame(table(subset(bp2,Abundance>=bcutoff)[,"id"])/length(unique(bp2$Sample)))
colnames(bp2)<-c("id","Prevalence")
bp2$"Group"<-"Bacteria"
prevalence$"Group"<-"Virus"

total.prevalence<-rbind(bp2,prevalence)

abb<-subset(All_abunt,Abundance>=vcutoff)[,c(1,2,3)]
VC_abundance_cast<-spread(abb,`sample id`,Abundance)
VC_abundance_cast[is.na(VC_abundance_cast)]<-0
rownames(VC_abundance_cast)<-VC_abundance_cast[,1]
VC_abundance_cast<-VC_abundance_cast[,-1]
# VC_abundance_cast<-VC_abundance_cast[which(rowSums(VC_abundance_cast) > (vcutoff*10)),]
# VC_abundance_cast<-t(VC_abundance_cast)
# VC_abundance_cast<-VC_abundance_cast[which(rowSums(VC_abundance_cast) > (vcutoff*10)),]
# VC_abundance_cast<-t(VC_abundance_cast)
VC_abundance_cast<-sweep(VC_abundance_cast,2,colSums(VC_abundance_cast),`/`)
# VC_abundance_cast<-t(VC_abundance_cast)
VC_shnnon<-data.frame(diversity(t(VC_abundance_cast), "shannon"))
bacprevalence2<-read_csv("/mnt/raid8/sunchuqing/2020_Human_phage/0001_BacMetagenomic/06_Abundance/genus_abundance_profile.csv")
bacprevalence2<-data.frame(bacprevalence2)
bp2<-bacprevalence2[,-1]
rownames(bp2)<-bacprevalence2$taxonomy
# bp2<-bp2[which(rowSums(bp2) > (bcutoff*10)),]
# bp2<-t(bp2)
# bp2<-bp2[which(rowSums(bp2) > (bcutoff*10)),]
# bp2<-t(bp2)
bp2<-bp2 %>%t() %>%melt() 
names(bp2)<-c("Sample","id","Abundance")
bp2$Abundance <- bp2$Abundance/100

bp2<-data.frame(subset(bp2,Abundance>=bcutoff))
bp2<-dcast(bp2,id~Sample)
bp2[is.na(bp2)]<-0
bac_shnnon<-data.frame(diversity(t(bp2[,-1]), "shannon"))
bac_shnnon$"Group"<-"Bac"

colnames(bac_shnnon)<-c("Shannon","Group")

VC_shnnon$"Group"<-"VC"

names(VC_shnnon)<-c("Shannon","Group")
# colors<-c("#ece2b9","#b5bbc3")
colors<-c("#e8d06e","#89a4c7")
p3<-rbind(bac_shnnon[,c(1,2)],VC_shnnon) %>%
  ggplot(aes(x=Group,y=Shannon))+
  geom_boxplot(aes(x=Group,y=Shannon,fill=Group),width=0.5)+
  # scale_fill_aaas()+
  scale_fill_manual(values =colors)+
  geom_signif(comparisons = list(c("Bac","VC")),step_increase = 0.1,
              map_signif_level = T,test = wilcox.test)+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  guides(fill=F,color=F)

cps2<-list(c("Bacteria",
             "Virus"))
p2<-total.prevalence %>% 
  mutate(id=as.character(id)) %>% 
  subset(Prevalence>0) %>% 
  ggplot(aes(x=Group,y=Prevalence))+
  geom_violin(aes(x=Group,y=Prevalence,fill=Group))+
  geom_boxplot(aes(x=Group,y=Prevalence),width=0.05,outlier.size = 0.4)+
  geom_signif(comparisons = cps2,step_increase = 0.1,
              map_signif_level = T,test = t.test)+
  # scale_fill_aaas()+
  scale_fill_manual(values =colors)+
  scale_color_manual(values =colors)+
  # scale_color_aaas()+
  scale_y_continuous(labels = scales::percent)+
  
  theme_bw()+
  theme(axis.text.x = element_blank())+
  labs(y="Prevalence",x="Group",title="",color="")
total.prevalence.1<-total.prevalence
A1<-(p2+p3) + plot_layout(guides = 'collect',tag_level = 'keep')

abb<-subset(All_abunt_sp,Abundance>=vcutoff)[,c(1,2,3)]
SP_abundance_cast<-spread(abb,`sample id`,Abundance)
SP_abundance_cast[is.na(SP_abundance_cast)]<-0
SP_abundance_cast<-data.frame(SP_abundance_cast)
rownames(SP_abundance_cast)<-SP_abundance_cast[,1]
SP_abundance_cast<-SP_abundance_cast[,-1]
# SP_abundance_cast<-SP_abundance_cast[which(rowSums(SP_abundance_cast) > (vcutoff*10)),]
# SP_abundance_cast<-t(SP_abundance_cast)
# SP_abundance_cast<-SP_abundance_cast[which(rowSums(SP_abundance_cast) > (vcutoff*10)),]
# SP_abundance_cast<-t(SP_abundance_cast)
SP_abundance_cast<-sweep(SP_abundance_cast,2,colSums(SP_abundance_cast),`/`)
SP_abundance_cast<-t(SP_abundance_cast) %>% melt()
names(SP_abundance_cast)<-c("Sample","id","Abundance")
abb<-SP_abundance_cast
prevalence<-data.frame(table(subset(All_abunt_sp,Abundance>=vcutoff)[,"id"])/length(unique(All_abunt_sp$`sample id`)))
colnames(prevalence)<-c("id","Prevalence")


abb<-subset(All_abunt_sp,Abundance>=vcutoff)[,c(1,2,3)]
SP_abundance_cast<-spread(abb,`sample id`,Abundance)
SP_abundance_cast[is.na(SP_abundance_cast)]<-0
SP_abundance_cast<-data.frame(SP_abundance_cast)
rownames(SP_abundance_cast)<-SP_abundance_cast[,1]
SP_abundance_cast<-SP_abundance_cast[,-1]
# SP_abundance_cast<-SP_abundance_cast[which(rowSums(SP_abundance_cast) > (vcutoff*10)),]
# SP_abundance_cast<-t(SP_abundance_cast)
# SP_abundance_cast<-SP_abundance_cast[which(rowSums(SP_abundance_cast) > (vcutoff*10)),]
# SP_abundance_cast<-t(SP_abundance_cast)
SP_abundance_cast<-sweep(SP_abundance_cast,2,colSums(SP_abundance_cast),`/`)
#SP_abundance_cast<-t(SP_abundance_cast)

SP_shnnon<-data.frame(diversity(t(SP_abundance_cast), "shannon"))
bacprevalence2<-read_csv("/mnt/raid8/sunchuqing/2020_Human_phage/0001_BacMetagenomic/06_Abundance/species_abundance_profile.csv")
bacprevalence2<-data.frame(bacprevalence2)
bp2<-bacprevalence2[,-1]
rownames(bp2)<-bacprevalence2$taxonomy
# bp2<-bp2[which(rowSums(bp2) > (bcutoff*10)),]
# bp2<-t(bp2)
# bp2<-bp2[which(rowSums(bp2) > (bcutoff*10)),]
# bp2<-t(bp2)
bp2<-bp2 %>%t() %>%melt() 
names(bp2)<-c("Sample","id","Abundance")
bp2$Abundance <- bp2$Abundance/100

bp2<-data.frame(subset(bp2,Abundance>=bcutoff))
bp2<-dcast(bp2,id~Sample)
bp2[is.na(bp2)]<-0
bac_shnnon<-data.frame(diversity(t(bp2[,-1]), "shannon"))
bac_shnnon$"Group"<-"Bac"

colnames(bac_shnnon)<-c("Shannon","Group")

SP_shnnon$"Group"<-"Vpops"
names(SP_shnnon)<-c("Shannon","Group")
colors<-c("#e8d06e","#89a4c7")

p1<-rbind(bac_shnnon,SP_shnnon) %>%
  ggplot(aes(x=Group,y=Shannon))+
  geom_boxplot(aes(x=Group,y=Shannon,fill=Group),width=0.5)+
  scale_fill_manual(values =colors)+
  geom_signif(comparisons = list(c("Bac","Vpops")),step_increase = 0.1,
              map_signif_level = T,test = t.test)+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  guides(fill=F,color=F)
# cbind(bac_shnnon,SP_shnnon)
bacprevalence2<-read_csv("/mnt/raid8/sunchuqing/2020_Human_phage/0001_BacMetagenomic/06_Abundance/species_abundance_profile.csv")
bacprevalence2<-data.frame(bacprevalence2)
bp2<-bacprevalence2[,-1]
rownames(bp2)<-bacprevalence2$taxonomy
# bp2<-bp2[which(rowSums(bp2) > (bcutoff*10)),]
# bp2<-t(bp2)
# bp2<-bp2[which(rowSums(bp2) > (bcutoff*10)),]
# bp2<-t(bp2)
bp2<-bp2 %>%t() %>%melt() 
names(bp2)<-c("Sample","id","Abundance")
bp2<-data.frame(table(subset(bp2,Abundance>=(bcutoff))[,"id"])/length(unique(bp2$Sample)))
colnames(bp2)<-c("id","Prevalence")
bp2$"Group"<-"Bacteria"
prevalence$"Group"<-"Virus"

total.prevalence<-data.frame(rbind(bp2,prevalence))


cp_sp<-list(c("Bacteria",
              "Virus"))
cp_sp<-list(unique(total.prevalence$Group))
p2<-total.prevalence %>% 
  mutate(id=as.character(id)) %>% 
  subset(Prevalence>0) %>% 
  ggplot(aes(x=Group,y=Prevalence))+
  geom_violin(aes(x=Group,y=Prevalence,fill=Group),adjust=2)+
  geom_boxplot(aes(x=Group,y=Prevalence),width=0.05,outlier.size = 0.4)+
  geom_signif(comparisons = cp_sp,step_increase = 0.1,
              map_signif_level = T,test = t.test)+
  scale_y_continuous(labels = scales::percent)+
  
  scale_fill_manual(values =colors)+
  scale_color_manual(values =colors)+
  theme_bw()+
  theme(axis.text.x = element_blank())+
  labs(y="Prevalence",x="Group",title="",color="")
A2<-(p2+p1) + plot_layout(guides = 'collect',tag_level = 'keep')
A<-(A2/A1)+ plot_layout(tag_level = 'new')
A

subset(total.prevalence.1,Group=="Bacteria")$Prevalence %>% median
subset(total.prevalence.1,Group=="Virus")$Prevalence %>% median
subset(total.prevalence,Group=="Bacteria")$Prevalence %>% median
subset(total.prevalence,Group=="Virus")$Prevalence %>% median

subset(total.prevalence.1,Group=="Bacteria")$Prevalence %>% mean
subset(total.prevalence.1,Group=="Virus")$Prevalence %>% mean
subset(total.prevalence,Group=="Bacteria")$Prevalence %>% mean
subset(total.prevalence,Group=="Virus")$Prevalence %>% mean

bac_shnnon$Shannon %>% median
VC_shnnon$Shannon %>% median
bac_shnnon$Shannon %>% median
SP_shnnon$Shannon %>% median


# 8. coverage -------------------------------------------------------------

annotated<-read.table("/mnt/raid8/sunchuqing/2022_OtherAna/Filtered.CHGV/HGV_annatation.tsv",sep = "\t",header=T)
annotated$"Group"<-"Novel"
annotated$"Group"[which(annotated$coverage>=0.7)]<-"Partially"
annotated$"Group"[which(annotated$coverage>=0.95)]<-"Identical"
annotated <- subset(annotated,contig %in% kept)
annotated %>%  
  ggplot(aes(x="Assembly",y=..count..,fill=Group))+
  geom_bar(width=0.7,position = "fill",color="white")+
  theme_void()+
  # theme(axis.x.test = element_blank()) +
  labs(y="",fill="")+
  # facet_wrap(~Assembly,ncol=2)+
  coord_polar(theta = 'y')+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_simpsons()


# 9. Top 10  --------------------------------------------------------------
cutoff <- 1e-6
prevalence<-data.frame(table(subset(All_abunt_sp,Abundance>=cutoff)[,"id"])/length(unique(All_abunt_sp$`sample id`)))
colnames(prevalence)<-c("id","Prevalence")
Top10 <- subset(VC,VC %in% paste("VC_",1:10,sep=""))
Top10 <- merge(Top10,prevalence[,c(1,2)],by.x="Contig",by.y="id",all.x=T)
Top10$Prevalence[is.na(Top10$Prevalence)] <- 0
Top10$Type[which(Top10$VC=="VC_1")] <- "VC_1"
Top10$Type[which(Top10$VC=="VC_4")] <- "VC_4"

Top10$Type[which(Top10$VC %in% paste("VC_",c(2,5,8,9,10),sep=""))] <- "crAssphage"
Top10$Type[which(Top10$VC %in% paste("VC_",c(3,6,7),sep=""))] <- "Gubaphage"
p1 <- Top10[,c("Contig","Type","VC","Length")] %>% 
  unique() %>% 
  mutate(VC=factor(VC,levels = paste("VC_",1:10,sep=""))) %>% 
  ggplot(aes(x=VC,color=Type,y=Length))+
  geom_boxplot()+
  scale_color_jco()+
  scale_y_log10()+
  theme_bw()+
  labs(x="")+
  theme(axis.text.x =element_blank())
p2 <- Top10[,c("Contig","Type","VC","Prevalence")] %>% 
  unique() %>%
  mutate(VC=factor(VC,levels = paste("VC_",1:10,sep=""))) %>% 
  ggplot(aes(x=VC,color=Type,y=Prevalence))+
  geom_boxplot()+
  scale_color_jco()+
  labs(x="")+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  theme(axis.text.x =element_blank())
Top10 <- merge(Top10,All_abunt_sp,by.x="Contig",by.y="id",all.x=T)
Top10$Abundance[is.na(Top10$Abundance)] <- 0
Top10.agg <- aggregate(Top10$Abundance,by=list(Top10$Contig,Top10$VC,Top10$Type),FUN=median)
colnames(Top10.agg) <- c("Contig","VC","Type","Median Abundance")
p3 <- Top10.agg %>% 
  mutate(VC=factor(VC,levels = paste("VC_",1:10,sep=""))) %>% 
  ggplot(aes(x=VC,color=Type,y=`Median Abundance`))+
  geom_boxplot()+
  scale_color_jco()+
  scale_y_log10()+
  labs(x="")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
(p1/p2/p3)+ plot_layout(guides = 'collect')

# 10. Enterotype ----------------------------------------------------------

library(reshape2)
library(readr)
library(ade4)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(vegan)#用于计算距离
library(edgeR)
abb<-subset(All_abunt,Abundance>=1e-5)[,c(1,2,3)]
VC_abundance_cast<-spread(abb,`sample id`,Abundance)
VC_abundance_cast[is.na(VC_abundance_cast)]<-0
rownames(VC_abundance_cast)<-VC_abundance_cast[,1]
VC_abundance_cast<-VC_abundance_cast[,-1]
VC_abundance_cast<-VC_abundance_cast[which(rowSums(VC_abundance_cast) > 0.001),]
VC_abundance_cast<-t(VC_abundance_cast)
VC_abundance_cast<-VC_abundance_cast[which(rowSums(VC_abundance_cast) > 0.001),]
VC_abundance_cast<-t(VC_abundance_cast)
VC_abundance_cast<-sweep(VC_abundance_cast,2,colSums(VC_abundance_cast),`/`)
data<-VC_abundance_cast 
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
KLD <- function(x,y) sum(x * log(x/y))
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}
# data.dist=dist.JSD(data)
# data.dist=vegdist(t(data),"euclidean")
data.dist=vegdist(t(data),"cao")
library(cluster)
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}
require(clusterSim)
nclusters=NULL

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
data.cluster=pam.clustering(data.dist, k=2)
data.cluster %>% table
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
obs.silhouette
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}
data.denoized=noise.removal(data, percent=0.01)
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F, col=c(6,5,4,2))
# s.class(obs.bet$ls, fac=as.factor(data.cluster), cstar=0, col=c(6,5,4,2))

obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, col=c(6,5,4,2))
#s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(6,5,4,2))

a<-subset(Abunt_top10,id=="Others")
label<-a[order(a$Abundance,decreasing = F),]$`sample id`
Abunt_top10$`sample id`<- factor(Abunt_top10$`sample id`, levels = label)

cluster.info<-cbind(colnames(data),data.cluster)
cluster.info<-data.frame(cluster.info)
names(cluster.info)<-c("sample id","Cluster")
abund.plot<-merge(Abunt_top10,cluster.info,by="sample id")
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colors<-c("Grey",getPalette(length(unique(Abunt_top10$id))))
aa <- subset(abund.plot,id=="Others")
label<-aa %>% arrange(-desc(Cluster),-desc(Abundance)) %>% dplyr::select("sample id") 
abund.plot$`sample id`<- factor(abund.plot$`sample id`, levels = as.character(label$`sample id`))
C1<-ggplot(abund.plot,aes(x=`sample id`,y=Abundance,fill=id,alluvium = id))+
  geom_bar(position ="fill",stat="identity",width = 0.5)+
  theme_bw()+
  guides(fill="none")+
  scale_fill_manual(values=colors)+
  #geom_alluvium(aes(fill = id))+
  facet_grid(~Cluster,scales = 'free_x',space="free_x")+
  geom_alluvium(aes(fill = id))+
  theme(axis.text.x = element_blank())+
  labs(title = "",y = "VC abundance",x="Sample ID")
C<-(C1+plot_spacer())+plot_layout(width=c(10,5))
C2<-(plot_spacer()+C1)+plot_layout(width=c(5,10))

(((B|A)+plot_layout(width=c(5,10)))/C)+
  plot_annotation(tag_levels = c('A','1'))+
  plot_layout(height=c(10,5))
ggsave(plot=C2,"/mnt/raid8/sunchuqing/R_proj/2020_Human_phage/Main.supp.plot/JustPlot/Main/Fig4C.pdf",width=15,height=5)
require(cowplot)
upper <- plot_grid(B, A, labels = c('A', 'B'), ncol = 2, rel_widths = c(5, 10))
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=T, col=c(6,5,4,2))
bottom_row <- C2
# then combine with the top row for final plot
fig4<-plot_grid(upper, bottom_row, labels = c('', 'C'), ncol = 1, rel_heights = c(10, 5),label_size = 14)




