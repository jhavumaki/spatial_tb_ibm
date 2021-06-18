rm(list=ls())
#setwd("~/Desktop/Postdoc/spatial_abm/IBM_hh/")
library(readr)
library(sp)
#library(GADMTools)
library(raster)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(pscl)
library(MASS)
library(boot)
library(MASS)
library(ks)
library(seewave)
library(spdep)
library(plot3Drgl)
#library(DescTools)
library(gridExtra)
library(ape)

dir<- "/gpfs/ysm/home/jh2822/project/spatial_abm/results_revisionsFINAL/"
dir2<-"/gpfs/ysm/home/jh2822/project/spatial_abm/hh_setup/"
ARRAYID=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(ARRAYID)

minLHS=1*(ARRAYID-1)+1
maxLHS=ARRAYID*1
numLHS=maxLHS-minLHS+1
study_time_frame=33 #33 months
intervention_time_frame=60

steady_state=11000
nb_depth=5
nb_depth2=3

#missing files
files<-list.files(dir) #list files in directory
LHS<-read.csv("lhs_10K_210419.csv")
rownames(LHS)<-LHS[,1]
LHS<-LHS[,-1]

# Neighbors list ----------------------------------------------------------
per <- getData('GADM', country='PER', level = 3)
smp=per[per$NAME_3=="San Martin de Porres",]
img<-raster("per_ppp_2009.tif")#("per_ppp_2009_UNadj.tif")
img2<-crop(img,extent(smp))
point_in_poly1<-mask(img2,smp)
point_in_polygon_df<-as(point_in_poly1, "SpatialPolygonsDataFrame")
w <- poly2nb(point_in_polygon_df, queen=T) #construct neighbors queen list, in the top corner, false for rook
#save(w, file="neighbor_list.Rdata")

#load("neighbor_list.Rdata", envir=.GlobalEnv)

#spatial scale; need to make sure these all allign

nb_list <- vector(mode = "list", length = length(w))
nb_list2 <- vector(mode = "list", length = length(w))
#this by default doesinclude own centroid 
for (kk in 1:length(w)){
  rr=1
  tempN<-w[[kk]]
  if (sum(tempN==0)==1){
    nb_list[[kk]]<-0
    
  } else {
while (rr<=nb_depth){
  nb_list[[kk]]<-unique(unlist(tempN, nb_list[[kk]]))
  tempN<-unlist(w[tempN])
  rr=rr+1
}
  }
}



for (kk in 1:length(w)){
  rr=1
  tempN<-w[[kk]]
  if (sum(tempN==0)==1){
    nb_list2[[kk]]<-0
    
  } else {
    while (rr<=nb_depth2){
      nb_list2[[kk]]<-unique(unlist(tempN, nb_list2[[kk]]))
      tempN<-unlist(w[tempN])
      rr=rr+1
    }
  }
}


# load files --------------------------------------------------------------
#model<-read.csv("outputHH7lhs53rng1862677429seed1none.csv")
data<-read.csv("case_centroid.csv")
location<-read.csv("centroid_location.csv")

location<-location[,-1]

data<-data[order(data$centroid),]
data<-data.frame(data %>% group_by(centroid) %>% tally())


strMatch<-"none"
strMatch2<-"output" #this corresponds to notifications
filesNo<-files[grepl(strMatch,files)]

filesNo1<-filesNo[grepl(strMatch2,filesNo)]
length(filesNo1)

slicer<-function(x) {
  df<-data.frame(model[,x:(x+study_time_frame-1)])
  colnames(df)<-steady_state+x:(x+study_time_frame-1)
  df<-rowSums(df)
return(df)
}
allKL<-data.frame( min_kl=numeric(), 
                   index=numeric(), 
                   fileID=numeric(), 
                   hh=numeric(), 
                   rng=numeric(), 
                   hh_atb=numeric(), 
                   detected=numeric(), 
                   moran_observed= numeric(),
                   moran_expected= numeric(),
                   moran_sd= numeric(),
                   moran_p= numeric(),
                   moran_observed_data= numeric(),
                   moran_expected_data= numeric(),
                   moran_sd_data= numeric(),
                   moran_p_data= numeric(),
                   moran_observed_mc= numeric(),
                   moran_p_mc= numeric(),
                   moran_observed_data_mc= numeric(),
                   moran_p_data_mc= numeric()
                   
                   
)


allobsMod5<-data.frame( obs_mod=numeric(), 
                        index=numeric(), 
                        fileID=numeric(), 
                        hh=numeric(), 
                        rng=numeric(), 
                        hh_atb=numeric(),
                        detected=numeric()
)

allobsMod3<-data.frame( obs_mod=numeric(), 
                        index=numeric(), 
                        fileID=numeric(), 
                        hh=numeric(), 
                        rng= numeric(), 
                        hh_atb=numeric(),
                        detected=numeric()
)





dataMat<-matrix(0, nrow=nrow(location), ncol= 5) 
dataMat[data$centroid,1]<-data$n
dataMat[,2]<-NA
dataMat[,3]<-1:nrow(location)
dataMat[,4]<-location[,1]
dataMat[,5]<-location[,2]
dataLocations<-data.frame(rep(dataMat[,4],dataMat[,1]), rep(dataMat[,5],dataMat[,1]))
#dataKDE<-  kde2d(x=dataLocations[,1], y=dataLocations[,2], n=500, 
#                 lims= c(min(location$longitude),max(location$longitude),min(location$latitude) , 
#                         max(location$latitude)))



dataKDE<- kde(dataLocations, eval.points=location[,1:2])
a<-data.frame(dataKDE$eval.points, abs(dataKDE$estimate))

png("calib_revisionsFINAL/DataKDE.png",  type="cairo")
ggplot(aes(a[,1], a[,2], color=a[,3]), data=a)+geom_point()+
  scale_color_gradientn(colours = terrain.colors(10))+
  xlab("Longitude")+ylab("Latitude") +
  theme(legend.title=element_blank())
dev.off()


for (i in minLHS:maxLHS) {
  strMatch3<-paste("lhs", i, "rng",sep="")
 
  filesNo2<-filesNo1[grepl(strMatch3,filesNo1)]
  print(filesNo2)
  if(length(filesNo2)>1){
    write.csv(filesNo2, paste("dup_file", i, ".csv",sep=""))
    filesNo2<-filesNo2[1]
  }  
  
   
  
  hh<-unlist(strsplit(filesNo2, "_none.csv",2))
  hh<-as.numeric(unlist(strsplit(hh,"HH"))[2])
  dist1<-read.csv(paste(dir2,"distance_table", hh, ".csv", sep="")) # this will be access via LHS
  dist1<-dist1[,-1]
  dist1[dist1$d==0,"d"]<-0.05
  
  rng<-unlist(strsplit(filesNo2, "rng"))[2]
  rng<-as.numeric(unlist(strsplit(rng, "seed"))[1])
 
  
  distTemp<-dist1
  distTemp<-distTemp[distTemp$d<=LHS[24,i],]
  distTemp$weights<-distTemp$d^(-LHS[28,i])
  
  
  
  model<-read.csv(paste(dir,filesNo2, sep="/")) # 11001 to 12000
  model=model[,-1]
  
  #model<-read.csv("/Users/joshuahavumaki/Downloads/output_state_nonelhs528rng555057635seedHH10_none.csv")
  
  resultsMat<-matrix(0, nrow=nrow(model), ncol= 1) 

  maxCol<-ncol(model)-(study_time_frame+intervention_time_frame-1)
  
  cols<-1:maxCol
  list<-lapply(cols, function(x) slicer(x)) #outputs a list with each element it's own 33 month slice
  
   kls<-matrix(0, nrow=length(list), ncol=3)
  obsMod5<-matrix(0, nrow=length(list), ncol=1)
  obsMod3<-matrix(0, nrow=length(list), ncol=1)
  

  
  for (j in 1:length(list)){
  

  dataMat[,2]<-list[[j]]
  resultsdf<-data.frame(dataMat)
  colnames(resultsdf)<-c("data", "model", "centroid", "longitude", "latitude")
  modelLocations<-data.frame(rep(resultsdf$longitude,resultsdf$model), rep(resultsdf$latitude,resultsdf$model ))

  tempState_ts<-read.csv(paste(dir,"/", "state_ts",unlist(strsplit(filesNo2, "output", "[[",2))[2], sep=""))
  
  temphh_atb<-read.csv(paste(dir,"/", "hh_atb",unlist(strsplit(filesNo2, "output_state", "[[",2))[2], sep=""))
  temphh_atb<-temphh_atb[,-1]
  
  temp_detected<-read.csv(paste(dir,"/", "detected",unlist(strsplit(filesNo2, "output_state_none", "[[",2))[2], sep=""))
  temp_detected<-temp_detected[,-1]
  
  #if (j==1){
  #save.image("temp.RData")
  #}
  #load("/Users/joshuahavumaki/Downloads/temp.RData")
  
  #temp_detected<-read.csv("/Users/joshuahavumaki/Downloads/detectedlhs880rng999873980seedHH2_none.csv")
  #temp_detected<-temp_detected[,-1]
  
  LTBI<-mean(tempState_ts$L[(steady_state+j):(steady_state+j+study_time_frame-1)])+ mean(tempState_ts$E[(steady_state+j):(steady_state+j+study_time_frame-1)])+ mean(tempState_ts$R[(steady_state+j):(steady_state+j+study_time_frame-1)])
  DetectedProportion<-sum(temp_detected[j:(j+study_time_frame-1),3])/sum(temp_detected[j:(j+study_time_frame-1),2])
    #(tail(tempState_ts$L,1) + tail(tempState_ts$E,1)+ tail(tempState_ts$R,1))
  
  if (!is.finite(LTBI) | !is.finite(DetectedProportion)) {
    
    kls[j,]<- 999
    
    obsMod5[j]<- 999
    obsMod3[j]<- 999
    } else if (nrow(modelLocations)>(0.25*nrow(dataLocations)+ nrow(dataLocations)) |
      nrow(modelLocations)<(nrow(dataLocations)-0.25*nrow(dataLocations) ) | LTBI<100000 | LTBI > 200000 | 
      DetectedProportion<0.8 |DetectedProportion>0.9 ){
    kls[j,]<- 999
   
    obsMod5[j]<- 999
    obsMod3[j]<- 999
  } else {
  #modelKDE<- kde2d(x=modelLocations[,1], y=modelLocations[,2], n=500, 
  #           lims= c(min(location$longitude),max(location$longitude),min(location$latitude) , 
  #           max(location$latitude)))
    
  
  modelKDE<- kde(modelLocations, eval.points=location[,1:2])
    
  #kls[j,]<-unlist(kl.dist(modelKDE$z,dataKDE$z,  base = exp(1)))   #D2 is kl distance 1st item wrt 2nd item, 2nd is reference/nors so use d2
  
  kls[j,]<-unlist(kl.dist(abs(modelKDE$estimate),abs(dataKDE$estimate),  base = exp(1)))

  obsMod5[j]<-mean(unlist((lapply(1:length(w), function(x) abs(sum(resultsdf$data[nb_list[[x]]])- sum(resultsdf$model[nb_list[[x]]]))))))
  
  obsMod3[j]<-mean(unlist((lapply(1:length(w), function(x) abs(sum(resultsdf$data[nb_list2[[x]]])- sum(resultsdf$model[nb_list2[[x]]]))))))
  
      
 
  
  }
  
  }
  


  # KL ----------------------------------------------------------------------
  klfinite<-kls[is.finite(kls[,2]),]
  
  minKLindx<-which(kls[,2]==min(klfinite[klfinite[,2]>0,2]))
  
  if (length(minKLindx)>1){
    indx<- sample(1:length(minKLindx),1)
    minKLindx<-minKLindx[indx]
  }
  allKL[i-minLHS+1,1]<-kls[minKLindx,2]
  allKL[i-minLHS+1,2]<-minKLindx
  allKL[i-minLHS+1,3]<-i
  allKL[i-minLHS+1,4]<-hh
  allKL[i-minLHS+1,5]<-rng
  allKL[i-minLHS+1,6]<-sum(temphh_atb[minKLindx:(minKLindx+study_time_frame-1),1])/sum(temphh_atb[minKLindx:(minKLindx+study_time_frame-1),2])
  allKL[i-minLHS+1,7]<-sum(temp_detected[minKLindx:(minKLindx+study_time_frame-1),3])/sum(temp_detected[minKLindx:(minKLindx+study_time_frame-1),2])
  
                            
  
  jj=minKLindx
  resultsdf$model<-NA
  resultsdf$model<-list[[jj]]
  modelLocationsKL<-data.frame(rep(resultsdf$longitude,resultsdf$model), rep(resultsdf$latitude,resultsdf$model ))
write.csv(modelLocationsKL, paste("calib_revisionsFINAL/modelKL_lhs",i,".csv", sep="" ))

wts<-matrix(0,nrow=length(w), ncol=length(w))
wts[cbind(distTemp$a,distTemp$b)]<-distTemp$weights

#save.image("bugTemp.RData")
#load("/Users/joshuahavumaki/Downloads/bugTemp.RData")

if(sum(resultsdf$model)>0 & sum(sum(wts))>0){
moranModel<-Moran.I(resultsdf$model, wts, alternative = "two.sided")
allKL[i-minLHS+1,"moran_observed"]<-moranModel$observed
allKL[i-minLHS+1,"moran_expected"]<-moranModel$expected
allKL[i-minLHS+1,"moran_sd"]<-moranModel$sd
allKL[i-minLHS+1,"moran_p"]<-moranModel$p.value

} else if (sum(resultsdf$model)==0 | sum(sum(wts))==0){
 
  allKL[i-minLHS+1,"moran_observed"]<-0
  allKL[i-minLHS+1,"moran_expected"]<-0
  allKL[i-minLHS+1,"moran_sd"]<-0
  allKL[i-minLHS+1,"moran_p"]<-0
   
}

#write.csv(resultsdf, "databug.csv")
#write.csv(wts, "wts.csv")
#resultsdf<-read.csv("/Users/joshuahavumaki/Downloads/databug.csv")
#wts<-read.csv("/Users/joshuahavumaki/Downloads/wts.csv")
#wts<-wts[,-1]
#save.image("bugTempWorking.RData")
#load("/Users/joshuahavumaki/Downloads/bugTempWorking.RData")
if (sum(sum(wts))>0){
moranData<-Moran.I(resultsdf$data, wts, alternative = "two.sided")

allKL[i-minLHS+1,"moran_observed_data"]<-moranData$observed
allKL[i-minLHS+1,"moran_expected_data"]<-moranData$expected
allKL[i-minLHS+1,"moran_sd_data"]<-moranData$sd
allKL[i-minLHS+1,"moran_p_data"]<-moranData$p.value

#moran MCMC
wts_listw<-mat2listw(wts)

moranModel<-moran.mc(resultsdf$model, wts_listw, nsim=999, alternative = "greater")
allKL[i-minLHS+1,"moran_observed_mc"]<-moranModel$statistic
allKL[i-minLHS+1,"moran_p_mc"]<-moranModel$p.value

moranData<-moran.mc(resultsdf$data, wts_listw, nsim=999, alternative = "greater")
allKL[i-minLHS+1,"moran_observed_data_mc"]<-moranData$statistic
allKL[i-minLHS+1,"moran_p_data_mc"]<-moranData$p.value
}  else if (sum(sum(wts))==0){ #this occurs when nothing is within spatial range
  
  allKL[i-minLHS+1,"moran_observed_data"]<-999
  allKL[i-minLHS+1,"moran_expected_data"]<-999
  allKL[i-minLHS+1,"moran_sd_data"]<-999
  allKL[i-minLHS+1,"moran_p_data"]<-999
  allKL[i-minLHS+1,"moran_observed_mc"]<-999
  allKL[i-minLHS+1,"moran_p_mc"]<-999
  allKL[i-minLHS+1,"moran_observed_data_mc"]<-999
  allKL[i-minLHS+1,"moran_p_data_mc"]<-999
  
}

  
  



#permutation
nsims=10001
allPerm_kl<-matrix(NA, nrow=10000)
for (simNr in 1:nsims){
  allPerm_kl[1]<-mean(unlist((lapply(1:length(w), function(x) abs(sum(resultsdf$data[nb_list[[x]]])- sum(resultsdf$model[nb_list[[x]]]))))))
  Model_Perm<-resultsdf$model[sample(1:length(w),length(w), replace=F)]
  allPerm_kl[simNr+1]<-mean(unlist((lapply(1:length(w), function(x) abs(sum(resultsdf$data[nb_list[[x]]])- sum(Model_Perm[nb_list[[x]]]))))))
}

write.csv(allPerm_kl, paste("calib_revisionsFINAL/allPerm_kl", i,".csv", sep="" ))



# ModObs5 ------------------------------------------------------------------
  min_obsMod5<-which(obsMod5==min(obsMod5))
  if (length(min_obsMod5)>1){
    indx<- sample(1:length(min_obsMod5),1)
    min_obsMod5<-min_obsMod5[indx]
  }
  allobsMod5[i-minLHS+1,1]<-obsMod5[min_obsMod5]
  allobsMod5[i-minLHS+1,2]<-min_obsMod5
  allobsMod5[i-minLHS+1,3]<-i
  allobsMod5[i-minLHS+1,4]<-hh
  allobsMod5[i-minLHS+1,5]<-rng
  allobsMod5[i-minLHS+1,6]<-sum(temphh_atb[min_obsMod5:(min_obsMod5+study_time_frame-1),1])/sum(temphh_atb[min_obsMod5:(min_obsMod5+study_time_frame-1),2])
  allobsMod5[i-minLHS+1,7]<-sum(temp_detected[min_obsMod5:(min_obsMod5+study_time_frame-1),3])/sum(temp_detected[min_obsMod5:(min_obsMod5+study_time_frame-1),2])
                            
  jj=min_obsMod5
  resultsdf$model<-NA
  resultsdf$model<-list[[jj]]
  modelLocationsMO5<-data.frame(rep(resultsdf$longitude,resultsdf$model), rep(resultsdf$latitude,resultsdf$model ))

  write.csv(modelLocationsMO5, paste("calib_revisionsFINAL/modelmo5_lhs", i,".csv", sep="" ))
  
  #permutation
  nsims=10001
  allPerm_om5<-matrix(NA, nrow=10000)
  for (simNr in 1:nsims){
    allPerm_om5[1]<-allobsMod5$obs_mod
    Model_Perm<-resultsdf$model[sample(1:length(w),length(w), replace=F)]
    allPerm_om5[simNr+1]<-mean(unlist((lapply(1:length(w), function(x) abs(sum(resultsdf$data[nb_list[[x]]])- sum(Model_Perm[nb_list[[x]]]))))))
  }
  
  write.csv(allPerm_om5, paste("calib_revisionsFINAL/allPerm_om5", i,".csv", sep="" ))
  
  
  # ModObs3 ------------------------------------------------------------------
  min_obsMod3<-which(obsMod3==min(obsMod3))
  if (length(min_obsMod3)>1){
    indx<- sample(1:length(min_obsMod3),1)
    min_obsMod3<-min_obsMod3[indx]
  }
  print(3)
  allobsMod3[i-minLHS+1,1]<-obsMod3[min_obsMod3]
  allobsMod3[i-minLHS+1,2]<-min_obsMod3
  allobsMod3[i-minLHS+1,3]<-i
  allobsMod3[i-minLHS+1,4]<-hh
  allobsMod3[i-minLHS+1,5]<-rng
  allobsMod3[i-minLHS+1,6]<-sum(temphh_atb[min_obsMod3:(min_obsMod3+study_time_frame-1),1])/sum(temphh_atb[min_obsMod3:(min_obsMod3+study_time_frame-1),2])
  allobsMod3[i-minLHS+1,7]<-sum(temp_detected[min_obsMod3:(min_obsMod3+study_time_frame-1),3])/sum(temp_detected[min_obsMod3:(min_obsMod3+study_time_frame-1),2])
  
  
  jj=min_obsMod3
  resultsdf$model<-NA
  resultsdf$model<-list[[jj]]
  modelLocationsMO3<-data.frame(rep(resultsdf$longitude,resultsdf$model), rep(resultsdf$latitude,resultsdf$model ))
  
  write.csv(modelLocationsMO3, paste("calib_revisionsFINAL/modelmo3_lhs", i,".csv", sep="" ))
  
  
  
  if (obsMod5[min_obsMod5] !=999){
    tempState_ts<-read.csv(paste(dir,"/", "state_ts",unlist(strsplit(filesNo2, "output", "[[",2))[2], sep=""))
    #image(kde2d(x=modelLocationsMO5[,1], y=modelLocationsMO5[,2], n=500, 
    #            lims= c(min(location$longitude),max(location$longitude),min(location$latitude) , 
    #            max(location$latitude))), main=paste("OM KDE 5:", obsMod5[min_obsMod5], "j=",min_obsMod5, "cases",nrow(modelLocationsMO5)),
    #      xlab="Longitude", ylab="Latitude")
    #image(kde2d(x=modelLocationsMO3[,1], y=modelLocationsMO3[,2], n=500, 
    #            lims= c(min(location$longitude),max(location$longitude),min(location$latitude) , 
    #                    max(location$latitude))), main=paste("OM KDE 3:", obsMod5[min_obsMod3], "j=",min_obsMod3, "cases",nrow(modelLocationsMO3)),
    #      xlab="Longitude", ylab="Latitude")
    #image(kde2d(x=modelLocationsKL[,1], y=modelLocationsKL[,2], n=500, 
    #            lims= c(min(location$longitude),max(location$longitude),min(location$latitude) , 
    #            max(location$latitude))), main=paste("KL KDE:", kls[minKLindx,2],  "j=",minKLindx, "cases",nrow(modelLocationsKL)),
    #      xlab="Longitude", ylab="Latitude")
    
    
    MO5KDE<- kde(modelLocationsMO5, eval.points=location[,1:2])
    
    save(MO5KDE, file=paste("calib_revisionsFINAL/mo5", i, "j",min_obsMod5, ".RData", sep=""))
    
    a<-data.frame(MO5KDE$eval.points, MO5KDE$estimate)
    colnames(a)<-c("longitude", "latitude" ,"estimate")
   p1<-   ggplot(aes(x=longitude, y=latitude, color=estimate), data=a)+geom_point()+
     scale_color_gradientn(colours = terrain.colors(10))+ ggtitle(paste("OM KDE 5:", obsMod5[min_obsMod5], "j=",min_obsMod5, "cases",nrow(modelLocationsMO5)))+
     xlab("Longitude")+ylab("Latitude") +
     theme(legend.title=element_blank())
   print({ png(paste("calib_revisionsFINAL/p1",  "i", i,".png", sep=""  ),  type="cairo")
 p1
   })
    dev.off()

   
   MO3KDE<- kde(modelLocationsMO3, eval.points=location[,1:2])
   
   save(MO3KDE, file=paste("calib_revisionsFINAL/mo3", i, "j",min_obsMod3, ".RData", sep=""))
   
   a1<-data.frame(MO3KDE$eval.points, MO3KDE$estimate)
p2<-ggplot(aes(a1[,1], a1[,2], color=a1[,3]), data=a1)+geom_point()+
  scale_color_gradientn(colours = terrain.colors(10))+ggtitle(paste("OM KDE 3:", obsMod5[min_obsMod3], "j=",min_obsMod3, "cases",nrow(modelLocationsMO3)))+
  xlab("Longitude")+ylab("Latitude") +
  theme(legend.title=element_blank())
   
   print({ png(paste("calib_revisionsFINAL/p2",  "i", i,".png", sep=""  ),  type="cairo")
p2
   })
   dev.off()
  
   
   
png(paste("calib_revisionsFINAL/p4",  "i", i,".png", sep=""  ),  type="cairo")
plot(tempState_ts$L, ylab="LL",xlab="months",main=paste("Steady State Check"))

   dev.off()
   
   KLKDE<- kde(modelLocationsKL, eval.points=location[,1:2])
   
   save(KLKDE, file=paste("calib_revisionsFINAL/kl", i, "j",minKLindx, ".RData", sep=""))
   
   a2<-data.frame(KLKDE$eval.points, abs(KLKDE$estimate))
     p3<-   ggplot(aes(a2[,1], a2[,2], color=a2[,3]), data=a2)+geom_point()+
       scale_color_gradientn(colours = terrain.colors(10))+ggtitle(paste("KL KDE:", kls[minKLindx,2],  "j=",minKLindx, "cases",nrow(modelLocationsKL)))+
       xlab("Longitude")+ylab("Latitude") +
       theme(legend.title=element_blank())
   print({  png(paste("calib_revisionsFINAL/p3",  "i", i,".png", sep=""  ),  type="cairo")
p3
  })
   dev.off()
  

   
   
   }  
}

write.csv(allKL,paste("calib_revisionsFINAL/allKL", ARRAYID,".csv", sep=""))
write.csv(allobsMod5, paste("calib_revisionsFINAL/allobsMod5", ARRAYID,".csv", sep=""))
print(4)
write.csv(allobsMod3, paste("calib_revisionsFINAL/allobsMod3", ARRAYID,".csv", sep=""))
