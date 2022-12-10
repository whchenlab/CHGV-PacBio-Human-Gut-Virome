get_VCabundance<-function(rpkm,vc,id,Contig){
  colnames(rpkm)<-c("sample","contig","rpkm")
  colnames(vc)<-c("VC","contig")
  x<-merge(rpkm,vc,by="contig") %>% subset(contig %in% Contig)

  y<-aggregate(x$rpkm,by=list(x$VC),FUN=sum)
  y$"Abundance"<-y$x/sum(y$x)
  y$"sample id"<-id
  colnames(y)<-c("id","RPKM","Abundance","sample id")
  species<-y[c("sample id","id","Abundance")]
  return(species)
}
get_SPabundance<-function(rpkm,Contig){
  colnames(rpkm)<-c("sample","contig","rpkm")
  x<-subset(rpkm,contig %in% Contig)
  x$"Abundance"<-x$rpkm/sum(x$rpkm)
  colnames(x)<-c("sample id","id","rpkm","Abundance")
  species<-x[c("sample id","id","Abundance")]
  return(species)
}
get_topn<-function(species,id,n){
  species <- species[order(species$Abundance, decreasing = TRUE), ]
  species_top10 <- species[1:n,]
  others<-1-sum(species_top10$Abundance)
  ot<-c(id,"Others",others)
  xx<-rbind(species_top10,ot)
  xx$Abundance<-as.numeric(xx$Abundance)
  return(xx)
}
plot_VCabundance<-function(sam,num,VC,Contig){
  library(readr)
  All_abunt<-data.frame("sample id"="","id"="","Abundance"="")[-1,]
  Abunt_top10<-data.frame("sample id"="","id"="","Abundance"="")[-1,]
  for(i in sam){
    if(i %in% All_abunt$`sample id` ) {
      next
    }
    ff<-paste("/mnt/raid8/sunchuqing/2020_Human_phage/09_RPKM/",i,"_NGS.RPKM.csv",sep="")
    temp<-read_csv(ff,col_names = F)
    tt<-get_VCabundance(temp,VC,i,Contig)
    tt_top<-get_topn(tt,i,num)
    #tt_top[,1]<-id
    All_abunt<-rbind(All_abunt,tt)
    Abunt_top10<-rbind(Abunt_top10,tt_top)
  }
  return(list(All_abunt,Abunt_top10))
}
plot_SPabundance<-function(sam,num,Contig){
  library(readr)
  All_abuntsp<-data.frame("sample id"="","id"="","Abundance"="")[-1,]
  Abunt_top10sp<-data.frame("sample id"="","id"="","Abundance"="")[-1,]
  for(i in sam){
    if(i %in% All_abuntsp$`sample id` ) {
      next
    }
    ff<-paste("/mnt/raid8/sunchuqing/2020_Human_phage/09_RPKM/",i,"_NGS.RPKM.csv",sep="")
    temp<-read_csv(ff,col_names = F)
    tt<-get_SPabundance(temp,Contig)
    tt_top<-get_topn(tt,i,num)
    #tt_top[,1]<-id
    All_abuntsp<-rbind(All_abuntsp,tt)
    Abunt_top10sp<-rbind(Abunt_top10sp,tt_top)
  }
  return(list(All_abuntsp,Abunt_top10sp))
}


