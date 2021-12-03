setwd("~/Desktop/Research/Dissertation Code and Data")
library(adehabitatHR)
library(rgdal)
library(raster)
library(geosphere)

## Home Range size ----
data1<-read.csv("HomeRange2020Diss.csv")
  
data1.T<-subset(data1, treatment=="transplant")
data1.C<-subset(data1, treatment=="control")
data1.L<-subset(data1, treatment=="longitudinal")

lizard.sp <- SpatialPoints(data1[c("long", "lat")])
proj4string(lizard.sp) = CRS("+init=epsg:32612")
lizard.sp <- spTransform(lizard.sp,CRS("+init=epsg:32612"))
# Transplant
data1.sp.T <- data1.T
data1.sp.T <- SpatialPointsDataFrame(coords=(cbind(data1.sp.T$long, data1.sp.T$lat)), data = data1.T["id"], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
data1.sp.T <- spTransform(data1.sp.T, CRS("+init=epsg:32612"))
plot(data1.sp.T, col = data1.sp.T$id)
# Control
data1.sp.C <- data1.C
data1.sp.C <- SpatialPointsDataFrame(coords=(cbind(data1.sp.C$long, data1.sp.C$lat)), data = data1.C["id"], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
data1.sp.C <- spTransform(data1.sp.C, CRS("+init=epsg:32612"))
plot(data1.sp.C, col = data1.sp.C$id)
# Longitudinal
data1.sp.L <- data1.L
data1.sp.L <- SpatialPointsDataFrame(coords=(cbind(data1.sp.L$long, data1.sp.L$lat)), data = data1.L["id"], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
data1.sp.L <- spTransform(data1.sp.L, CRS("+init=epsg:32612"))
plot(data1.sp.L, col = data1.sp.L$id)

# 100%
# Transplant
lizard.mcp.100.T <- mcp(data1.sp.T[, "id"], percent = 100, unin = "m", unout = "m" )
as.data.frame(lizard.mcp.100.T)
plot(lizard.mcp.100.T)
plot(data1.sp.T)
# Control
lizard.mcp.100.C <- mcp(data1.sp.C[, "id"], percent = 100, unin = "m", unout = "m" )
as.data.frame(lizard.mcp.100.C)
plot(lizard.mcp.100.C)
plot(data1.sp.C)
# Longitudinal
lizard.mcp.100.L <- mcp(data1.sp.L[, "id"], percent = 100, unin = "m", unout = "m" )
as.data.frame(lizard.mcp.100.L)
plot(lizard.mcp.100.L)
plot(data1.sp.L)

lizard.mcp.100.L<-as.data.frame(lizard.mcp.100.L)
lizard.mcp.100.L$Treatment<-"Longitudinal"
lizard.mcp.100.T<-as.data.frame(lizard.mcp.100.T)
lizard.mcp.100.T$Treatment<-"Transplant"
lizard.mcp.100.C<-as.data.frame(lizard.mcp.100.C)
lizard.mcp.100.C$Treatment<-"Control"
lizard.mcp.100<-rbind(lizard.mcp.100.L, lizard.mcp.100.T, lizard.mcp.100.C)
colnames(lizard.mcp.100)[2]<-"HR100"

# 95%
# Transplant
lizard.mcp.95.T <- mcp(data1.sp.T[, "id"], percent = 95, unin = "m", unout = "m" )
as.data.frame(lizard.mcp.95.T)
plot(lizard.mcp.95.T)
plot(data1.sp.T)
# Control
lizard.mcp.95.C <- mcp(data1.sp.C[, "id"], percent = 95, unin = "m", unout = "m" )
as.data.frame(lizard.mcp.95.C)
plot(lizard.mcp.95.C)
plot(data1.sp.C)
# Longitudinal
lizard.mcp.95.L <- mcp(data1.sp.L[, "id"], percent = 95, unin = "m", unout = "m" )
as.data.frame(lizard.mcp.95.L)
plot(lizard.mcp.95.L)
plot(data1.sp.L)

lizard.mcp.95.L<-as.data.frame(lizard.mcp.95.L)
lizard.mcp.95.L$Treatment<-"Longitudinal"
lizard.mcp.95.T<-as.data.frame(lizard.mcp.95.T)
lizard.mcp.95.T$Treatment<-"Transplant"
lizard.mcp.95.C<-as.data.frame(lizard.mcp.95.C)
lizard.mcp.95.C$Treatment<-"Control"
lizard.mcp.95<-rbind(lizard.mcp.95.L, lizard.mcp.95.T, lizard.mcp.95.C)
colnames(lizard.mcp.95)[2]<-"HR95"

# Home Range 95% and 100% MCP
lizardHR<-merge(lizard.mcp.95, lizard.mcp.100)
#write.csv(lizardHR, "HomeRange.csv", row.names=F)
#writeOGR(lizard.mcp.95.L, dsn = ".", layer = "Longitudinal.95", driver="ESRI Shapefile")
# Maps were created by extracting shapefiles (Example in line 93) to use in QGIS

## Daily Distance Traveled  ----
data1$date<-as.Date(data1$date, format='%m/%d/%y')
data1<-data1[order(data1$date),]

p2002<-subset(data1, id=="2002")
p2002<-as.matrix(p2002[, c(4,3)])
a=as.data.frame(sum(distGeo(p2002), na.rm=T)) # change the function to mean() for each individual to calculate mean distance traveled
colnames(a)<-"SumDistance"
a$id<-"2002"

p2005<-subset(data1, id=="2005")
p2005<-as.matrix(p2005[, c(4,3)])
b=as.data.frame(sum(distGeo(p2005), na.rm=T))
colnames(b)<-"SumDistance"
b$id<-"2005"

p2008<-subset(data1, id=="2008")
p2008<-as.matrix(p2008[, c(4,3)])
c=as.data.frame(sum(distGeo(p2008), na.rm=T))
colnames(c)<-"SumDistance"
c$id<-"2008"

p2010<-subset(data1, id=="2010")
p2010<-as.matrix(p2010[, c(4,3)])
d=as.data.frame(sum(distGeo(p2010), na.rm=T))
colnames(d)<-"SumDistance"
d$id<-"2010"

p2011<-subset(data1, id=="2011")
p2011<-as.matrix(p2011[, c(4,3)])
e=as.data.frame(sum(distGeo(p2011), na.rm=T))
colnames(e)<-"SumDistance"
e$id<-"2011"

p2012<-subset(data1, id=="2012")
p2012<-as.matrix(p2012[, c(4,3)])
f=as.data.frame(sum(distGeo(p2012), na.rm=T))
colnames(f)<-"SumDistance"
f$id<-"2012"

p2013<-subset(data1, id=="2013")
p2013<-as.matrix(p2013[, c(4,3)])
g=as.data.frame(sum(distGeo(p2013), na.rm=T))
colnames(g)<-"SumDistance"
g$id<-"2013"

p2015<-subset(data1, id=="2015")
p2015<-as.matrix(p2015[, c(4,3)])
h=as.data.frame(sum(distGeo(p2015), na.rm=T))
colnames(h)<-"SumDistance"
h$id<-"2015"

p2016<-subset(data1, id=="2016")
p2016<-as.matrix(p2016[, c(4,3)])
i=as.data.frame(sum(distGeo(p2016), na.rm=T))
colnames(i)<-"SumDistance"
i$id<-"2016"

p2019<-subset(data1, id=="2019")
p2019<-as.matrix(p2019[, c(4,3)])
j=as.data.frame(sum(distGeo(p2019), na.rm=T))
colnames(j)<-"SumDistance"
j$id<-"2019"

p2020<-subset(data1, id=="2020")
p2020<-as.matrix(p2020[, c(4,3)])
k=as.data.frame(sum(distGeo(p2020), na.rm=T))
colnames(k)<-"SumDistance"
k$id<-"2020"

p2021<-subset(data1, id=="2021")
p2021<-as.matrix(p2021[, c(4,3)])
l=as.data.frame(sum(distGeo(p2021), na.rm=T))
colnames(l)<-"SumDistance"
l$id<-"2021"

p2022<-subset(data1, id=="2022")
p2022<-as.matrix(p2022[, c(4,3)])
m=as.data.frame(sum(distGeo(p2022), na.rm=T))
colnames(m)<-"SumDistance"
m$id<-"2022"

p2024<-subset(data1, id=="2024")
p2024<-as.matrix(p2024[, c(4,3)])
n=as.data.frame(sum(distGeo(p2024), na.rm=T))
colnames(n)<-"SumDistance"
n$id<-"2024"

p2028<-subset(data1, id=="2028")
p2028<-as.matrix(p2028[, c(4,3)])
o=as.data.frame(sum(distGeo(p2028), na.rm=T))
colnames(o)<-"SumDistance"
o$id<-"2028"

p2033<-subset(data1, id=="2033")
p2033<-as.matrix(p2033[, c(4,3)])
p=as.data.frame(sum(distGeo(p2033), na.rm=T))
colnames(p)<-"SumDistance"
p$id<-"2033"

p2034<-subset(data1, id=="2034")
p2034<-as.matrix(p2034[, c(4,3)])
q=as.data.frame(sum(distGeo(p2034), na.rm=T))
colnames(q)<-"SumDistance"
q$id<-"2034"

p2035<-subset(data1, id=="2035")
p2035<-as.matrix(p2035[, c(4,3)])
r=as.data.frame(sum(distGeo(p2035), na.rm=T))
colnames(r)<-"SumDistance"
r$id<-"2035"

TotalDistance<-rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r)



