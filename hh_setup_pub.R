rm(list=ls())
setwd("hh_setup")
library(raster)
library(GADMTools)
library(ggplot2)
library(sp)
library(triangle)
library(tiff)

# distribution of Peru HH sizes from http://data.un.org/Data.aspx?d=POP&f=tableCode:50------------------------------------------------
peru = c(
  "1" =794661,
  "2"=943300,
  "3"=1282876,
  "4"=1340476, 
  "5"=980549,
  "6+"=1412212
)

#empirical household size distribution
pdf("peru_hh_size.pdf", height=5)
barplot(peru, main="Peru Household distribution", xlab="Household Sizes", ylab="Count")
dev.off()


# Spatial area setup ------------------------------------------------------
per <- getData('GADM', country='PER', level = 3)
smp=per[per$NAME_3=="San Martin de Porres",]
class(smp)

img<-raster("per_ppp_2009.tif")#download from worldpop.org; unadjusted population density for 2009
#plot(img)

img2<-crop(img,extent(smp))
point_in_poly<-mask(img2,smp)
point_in_poly_df<-as(point_in_poly, "SpatialPixelsDataFrame")
point_in_poly_df <- as.data.frame(point_in_poly_df)
colnames(point_in_poly_df) <- c("value", "x", "y")

centroid <- unname(rasterToPoints(point_in_poly))
sum(centroid[,3])
mean(centroid[,3])

for (j in 1:10){
hh_data<-data.frame(
x=numeric(),
y=numeric(),
size=numeric(),
centroid =numeric()
)

#centroid[,3]<-centroid[,3]*4

#sample household sizes for Modeled population
for (i in 1:nrow(centroid)){
targetPop<-round(centroid[i,3])
HHs<-c()
while(sum(HHs)<targetPop){
#HHstemp<-as.integer(rtriangle(n=1, a = 1, b = 9, c = 4)) #http://data.un.org/Data.aspx?d=POP&f=tableCode:50

HHstemp<-sample(1:6, 1, prob=peru/sum(peru), replace=T)
if (HHstemp==6){
  
  HHstemp<-as.integer(rtriangle(n=1, a = 6, b = 10, c = 6))
}


if (sum(c(HHs, HHstemp))>targetPop){
  break 
}else if (sum(c(HHs, HHstemp))<=targetPop){
HHs<-c(HHs, HHstemp)
}
}

hh_temp<-data.frame(centroid[i,1],centroid[i,2],HHs, i)
colnames(hh_temp)<-c("x","y", "size", "centroid")
hh_data<-rbind(hh_data,hh_temp)
}
sum(HHs)

sum(hh_data$size)
sum(round(centroid[,3]))
a<-data.frame(table(hh_data$size))
a[,1]<-as.numeric(as.character(a[,1]))
a[,2]<-as.numeric(as.character(a[,2]))
a[6:9,1]<-6

model_hh_size<-c(a[1:5,2],  sum(a[6:nrow(a),2]))
model_hh_size/sum(model_hh_size)
peru/sum(peru)
a<-rep(a[,1], a[,2])
unique(a)
hist(a, breaks=length(unique(a))+1)
hh_data_cent<-aggregate(hh_data$size, by=list(Category=hh_data$centroid), FUN=sum)

a<-hh_data[match( hh_data_cent$Category, hh_data$centroid),c("x", "y")]
hh_data_cent<-data.frame(a,hh_data_cent )  
colnames(hh_data_cent)<-c("x","y", "centroid", "size")
ggplot(aes(x=x, y=y, fill=size), data=hh_data_cent)+geom_tile() +scale_fill_gradient(low="blue", high="red") 


hist(hh_data$size, breaks=10)

peru/sum(peru)

# distance calculation ----------------------------------------------------
mat<-matrix(0, ncol=nrow(centroid), nrow=nrow(centroid))
e<-upper.tri(mat, diag = FALSE) #to divide matrix
pair_df<-as.data.frame(which(e==1, arr.ind = T)) #to get HH connections
colnames(pair_df)<-c("a","b")
rm(e)
rm(mat)


long1 <- centroid[pair_df$a,1]*pi/180
lat1 <- centroid[pair_df$a,2]*pi/180  
long2 <- centroid[pair_df$b,1]*pi/180
lat2 <- centroid[pair_df$b,2]*pi/180  

dlong = (long2 - long1)
dlat  = (lat2 - lat1)
# Haversine formula:
R = 6371 #in km
a = sin(dlat/2)*sin(dlat/2) + cos(lat1)*cos(lat2)*sin(dlong/2)*sin(dlong/2)
c = 2 * atan2( sqrt(a), sqrt(1-a) )
pair_df$d = R * c 

sum(pair_df$a==2 & pair_df$b==1)
sum(pair_df$a==1 & pair_df$b==2)
revDf<-data.frame(pair_df$b,pair_df$a,  pair_df$d)  #not symetric so need to 1 to 2 and 2 to 1 
selfDF<-data.frame(unique(hh_data$centroid), unique(hh_data$centroid),0)
colnames(revDf)<-c("a", "b", "d")
colnames(selfDF)<-colnames(revDf)
pair_df<-rbind(pair_df,revDf,selfDF )
sum(pair_df$a==2 & pair_df$b==1)
sum(pair_df$a==1 & pair_df$b==2)
which(pair_df$a==2 & pair_df$b==1 )
which(pair_df$a==1 & pair_df$b==2 )

write.csv(pair_df,paste("distance_table", j,".csv", sep=""))
write.csv(hh_data, paste("hh_df", j,".csv", sep="") ) 

}
  

#modeled HH distribution
for (pp in 1:10){
hh_df<-read.csv(paste("/Volumes/Extreme SSD/spatial_abm/hh_setup/hh_df",pp, ".csv", sep=""))
pdf("model_hh_size.pdf", height=5)
HH_sizes<-ifelse(hh_df$size>5, "6+", hh_df$size)
HH_sizes<-table(HH_sizes)
barplot(HH_sizes, main="Modeled Household distribution", xlab="Household Sizes", ylab="Count")
dev.off()
}
