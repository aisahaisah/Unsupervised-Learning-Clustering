library(xlsx)
library(factoextra)
library(ggplot2)
library(clValid)
library(cluster)
library(FactoMineR)
library(psych)
library(fclust)
###Input Data###
dt<-read.xlsx(file.choose(),1)
###Memilih Data X###
dataX<-dt[,-1]
###Descriptive Statistics###
summary(dataX)

###Test Korelasi Multivariat##
mcor<-cor(dataX)
bartlet.test <- cortest.bartlett(mcor,n=154,diag=TRUE)

####Principal Component Analysis#####
x.pca<-PCA(dataX, scale.unit = TRUE, ncp = 8, ind.sup = NULL, quanti.sup = NULL, quali.sup = NULL, row.w = NULL, col.w = NULL, graph = TRUE, axes = c(1,2))
eigenvalue<-x.pca$eig
x.reduksi<-x.pca$ind$coord
xdim4<-x.reduksi[,-(5:8)]

###Preprocessing data######
###Normalisasi###
normalize <- function(M) {
  #center data
  means = apply(M,2,mean)
  Xnorm = t(apply(M,1,function(x) {x-means}))
  Xnorm
}
encode <- function(M) {
  # put on hypershpere
  mins = apply(M,2,min)
  maxs = apply(M,2,max)
  ranges = maxs-mins
  Xnorm = t(apply(M,1,function(x) { 2*(x-mins)/ranges-1}))
  Xnorm = t(apply(Xnorm, 1, function(x) { x/norm(x,type="2")}))
  Xnorm
}

###############4 DIMENSI#####################
dataX.4=xdim4
dataX.norm4=normalize(dataX.4)
dataX.norm4z=encode(dataX.norm4)

###K-Means Clustering###
fkmeans=function(data,k){
  kmeanss=kmeans(data,k,iter.max=200,nstart = 10)
  print(cbind.data.frame(Wilayah=dt$Nama.Wilayah,kluster=kmeanss$cluster))
}
fkmeans2<-fkmeans(dataX.norm4z,2)
fkmeans3<-fkmeans(dataX.norm4z,3)
fkmeans4<-fkmeans(dataX.norm4z,4)
fkmeans5<-fkmeans(dataX.norm4z,5)
fkmeans6<-fkmeans(dataX.norm4z,6)

##Plot KMeans###
kmeans2=kmeans(dataX.norm4z,2,iter.max=200,nstart = 10)
rownames(dataX.norm4z)=dt$Nama.Wilayah
kmeans3=kmeans(dataX.norm4z,3,iter.max=200,nstart = 10)
kmeans4=kmeans(dataX.norm4z,4,iter.max=200,nstart = 10)
kmeans5=kmeans(dataX.norm4z,5,iter.max=200,nstart = 10)
kmeans6=kmeans(dataX.norm4z,6,iter.max=200,nstart = 10)
fviz_cluster(kmeans2, data = dataX.norm4z,
             palette = c("#00FFFF","#2E9FDF", "#E7B800"),
             ggtheme = theme_minimal(),
             main = "Partitioning Clustering Plot")

plotk2=fviz_cluster(kmeans2, data = dataX.norm4z, choose.vars = NULL, stand = TRUE,
                    axes = c(1, 2), geom = c("point", "text"), repel = FALSE,
                    show.clust.cent = TRUE, ellipse = TRUE, ellipse.type = "convex",
                    ellipse.level = 0.95, ellipse.alpha = 0.2, shape = NULL,
                    pointsize = 1.5, labelsize = 8, main = "Cluster plot", xlab = NULL,
                    ylab = NULL, outlier.color = "black", outlier.shape = 19,
                    ggtheme = theme_grey())

plotk3=fviz_cluster(kmeans3, data = dataX.norm4z, choose.vars = NULL, stand = TRUE,
                    axes = c(1, 2), geom = c("point", "text"), repel = FALSE,
                    show.clust.cent = TRUE, ellipse = TRUE, ellipse.type = "convex",
                    ellipse.level = 0.95, ellipse.alpha = 0.2, shape = NULL,
                    pointsize = 1.5, labelsize = 8, main = "Cluster plot", xlab = NULL,
                    ylab = NULL, outlier.color = "black", outlier.shape = 19,
                    ggtheme = theme_grey())

plotk4=fviz_cluster(kmeans4, data = dataX.norm4z, choose.vars = NULL, stand = TRUE,
                    axes = c(1, 2), geom = c("point", "text"), repel = FALSE,
                    show.clust.cent = TRUE, ellipse = TRUE, ellipse.type = "convex",
                    ellipse.level = 0.95, ellipse.alpha = 0.2, shape = NULL,
                    pointsize = 1.5, labelsize = 8, main = "Cluster plot", xlab = NULL,
                    ylab = NULL, outlier.color = "black", outlier.shape = 19,
                    ggtheme = theme_grey())

plotk5=fviz_cluster(kmeans5, data = dataX.norm4z, choose.vars = NULL, stand = TRUE,
                    axes = c(1, 2), geom = c("point", "text"), repel = FALSE,
                    show.clust.cent = TRUE, ellipse = TRUE, ellipse.type = "convex",
                    ellipse.level = 0.95, ellipse.alpha = 0.2, shape = NULL,
                    pointsize = 1.5, labelsize = 8, main = "Cluster plot", xlab = NULL,
                    ylab = NULL, outlier.color = "black", outlier.shape = 19,
                    ggtheme = theme_grey())

plotk6=fviz_cluster(kmeans6, data = dataX.norm4z, choose.vars = NULL, stand = TRUE,
                    axes = c(1, 2), geom = c("point", "text"), repel = FALSE,
                    show.clust.cent = TRUE, ellipse = TRUE, ellipse.type = "convex",
                    ellipse.level = 0.95, ellipse.alpha = 0.2, shape = NULL,
                    pointsize = 1.5, labelsize = 8, main = "Cluster plot", xlab = NULL,
                    ylab = NULL, outlier.color = "black", outlier.shape = 19,
                    ggtheme = theme_grey())

###Mencari k optimal###
##Cara 1 Connectivity dan Silhouette###
fviz_nbclust(dataX.norm4z, kmeans, method = "silhouette")
perb.cluster <- clValid(dataX.norm4z, nClust = 2:6, 
                        clMethods = c("kmeans"), validation = "internal")
summary(perb.cluster)

###Cara 2 PseudoF###
icdrate = function(Data, nc, c)
{
  n = dim(Data)[1]  #jumlah baris
  p = dim(Data)[2]  #jumlah kolom
  X = Data[,1:(p-1)]
  Group = Data[,p]
  
  p = dim(X)[2]
  Mean.X = matrix(ncol = p, nrow = (nc+1))
  for (i in 1:nc)
  {
    for (j in 1:p)
    {
      Mean.X[i,j] = mean(X[which(Group==i),j])
      Mean.X[(nc+1),j] = mean(X[,j])
    }
  }
  
  SST = matrix(ncol=p, nrow=n)
  for (i in 1:n)
  {
    for (j in 1:p)
    {
      SST[i,j] = (X[i,j] - Mean.X[(nc+1),j])^2
    }
  }
  SST = sum(sum(SST))
  
  SSE = matrix(ncol=p, nrow=n)
  for (i in 1:n)
  {
    for (j in 1:p)
    {
      for (k in 1:nc)
      {
        if (Group[i]==k)
        {
          SSE[i,j] = (X[i,j] - Mean.X[k,j])^2
        }
      }
    }
  }
  SSE = sum(sum(SSE))
  
  Rsq = (SST-SSE)/SST
  pseudof = (Rsq/(c-1))/((1-Rsq)/(nc-c))
  icdrate = 1-Rsq
  list(Rsq = Rsq, pseudof = pseudof, icdrate = icdrate)
}
datainputps2=cbind(dataX.norm4z,fkmeans2$kluster)
datainputps3=cbind(dataX.norm4z,fkmeans3$kluster)
datainputps4=cbind(dataX.norm4z,fkmeans4$kluster)
datainputps5=cbind(dataX.norm4z,fkmeans5$kluster)
datainputps6=cbind(dataX.norm4z,fkmeans6$kluster)
ujips2 = icdrate(datainputps2,154,2)
ujips3 = icdrate(datainputps3,154,3)
ujips4 = icdrate(datainputps4,154,4)
ujips5 = icdrate(datainputps5,154,5)
ujips6 = icdrate(datainputps6,154,6)
kOptimKM=cbind.data.frame(k2=c(ujips2$pseudof,ujips2$icdrate),k3=c(ujips3$pseudof,ujips3$icdrate),k4=c(ujips3$pseudof,ujips4$icdrate),k5=c(ujips5$pseudof,ujips5$icdrate),k6=c(ujips6$pseudof,ujips6$icdrate)) 
rownames(kOptim)=c("pseudof","icdrate")

####GustafsonKessell###
library(ppclust)
# Initialize the prototype matrix using K-means++
v2 <- inaparc::kmpp(dataX.norm4z, k=2)$v
v3 <- inaparc::kmpp(dataX.norm4z, k=3)$v
v4 <- inaparc::kmpp(dataX.norm4z, k=4)$v
v5 <- inaparc::kmpp(dataX.norm4z, k=5)$v
v6 <- inaparc::kmpp(dataX.norm4z, k=6)$v
v7 <- inaparc::kmpp(dataX.norm4z, k=7)$v
v8 <- inaparc::kmpp(dataX.norm4z, k=8)$v
v9 <- inaparc::kmpp(dataX.norm4z, k=9)$v
v10 <- inaparc::kmpp(dataX.norm4z, k=10)$v
# Initialize the memberships degrees matrix
u2 <- inaparc::imembrand(nrow(dataX.norm4z), k=2)$u
u3 <- inaparc::imembrand(nrow(dataX.norm4z), k=3)$u
u4 <- inaparc::imembrand(nrow(dataX.norm4z), k=4)$u
u5 <- inaparc::imembrand(nrow(dataX.norm4z), k=5)$u
u6 <- inaparc::imembrand(nrow(dataX.norm4z), k=6)$u
u7 <- inaparc::imembrand(nrow(dataX.norm4z), k=7)$u
u8 <- inaparc::imembrand(nrow(dataX.norm4z), k=8)$u
u9 <- inaparc::imembrand(nrow(dataX.norm4z), k=9)$u
u10 <- inaparc::imembrand(nrow(dataX.norm4z), k=10)$u
##hasil clustering GK###
gkpfcm.result2=gkpfcm(dataX.norm4z, centers=v2, memberships=u2, m=2, eta=2, K=1, a=1, b=1, 
                      dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
                      alginitu="imembrand", nstart=1, iter.max=1000, con.val=1e-09, 
                      fixcent=FALSE, fixmemb=FALSE, stand=FALSE)
gkpfcm.result3=gkpfcm(dataX.norm4z, centers=v3, memberships=u3, m=2, eta=2, K=1, a=1, b=1, 
                      dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
                      alginitu="imembrand", nstart=1, iter.max=1000, con.val=1e-09, 
                      fixcent=FALSE, fixmemb=FALSE, stand=FALSE)
gkpfcm.result4=gkpfcm(dataX.norm4z, centers=v4, memberships=u4, m=2, eta=2, K=1, a=1, b=1, 
                      dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
                      alginitu="imembrand", nstart=1, iter.max=1000, con.val=1e-09, 
                      fixcent=FALSE, fixmemb=FALSE, stand=FALSE)
gkpfcm.result5=gkpfcm(dataX.norm4z, centers=v5, memberships=u5, m=2, eta=2, K=1, a=1, b=1, 
                      dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
                      alginitu="imembrand", nstart=1, iter.max=1000, con.val=1e-09, 
                      fixcent=FALSE, fixmemb=FALSE, stand=FALSE)
gkpfcm.result6=gkpfcm(dataX.norm4z, centers=v6, memberships=u6, m=2, eta=2, K=1, a=1, b=1, 
                      dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
                      alginitu="imembrand", nstart=1, iter.max=1000, con.val=1e-09, 
                      fixcent=FALSE, fixmemb=FALSE, stand=FALSE)
gkpfcm.result7=gkpfcm(dataX.norm4z, centers=v7, memberships=u7, m=2, eta=2, K=1, a=1, b=1, 
                      dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
                      alginitu="imembrand", nstart=1, iter.max=1000, con.val=1e-09, 
                      fixcent=FALSE, fixmemb=FALSE, stand=FALSE)
gkpfcm.result8=gkpfcm(dataX.norm4z, centers=v8, memberships=u8, m=2, eta=2, K=1, a=1, b=1, 
                      dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
                      alginitu="imembrand", nstart=1, iter.max=1000, con.val=1e-09, 
                      fixcent=FALSE, fixmemb=FALSE, stand=FALSE)
gkpfcm.result9=gkpfcm(dataX.norm4z, centers=v9, memberships=u9, m=2, eta=2, K=1, a=1, b=1, 
                      dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
                      alginitu="imembrand", nstart=1, iter.max=1000, con.val=1e-09, 
                      fixcent=FALSE, fixmemb=FALSE, stand=FALSE)
gkpfcm.result10=gkpfcm(dataX.norm4z, centers=v10, memberships=u10, m=2, eta=2, K=1, a=1, b=1, 
                      dmetric="sqeuclidean", pw = 2, alginitv="kmpp", 
                      alginitu="imembrand", nstart=1, iter.max=1000, con.val=1e-09, 
                      fixcent=FALSE, fixmemb=FALSE, stand=FALSE)
mem.gk2<-cbind.data.frame(dt$Nama.Wilayah,gkpfcm.result2$cluster)
mem.gk3<-cbind.data.frame(dt$Nama.Wilayah,gkpfcm.result3$cluster)
mem.gk4<-cbind.data.frame(dt$Nama.Wilayah,gkpfcm.result4$cluster)
mem.gk5<-cbind.data.frame(dt$Nama.Wilayah,gkpfcm.result5$cluster)
mem.gk6<-cbind.data.frame(dt$Nama.Wilayah,gkpfcm.result6$cluster)

membershipGK<-cbind(dt$Nama.Wilayah,gkpfcm.result2$cluster,gkpfcm.result3$cluster,gkpfcm.result4$cluster,gkpfcm.result5$cluster,gkpfcm.result6$cluster)

###Silhouette###
silgk2 <- SIL.F(dataX.norm4z,gkpfcm.result2$u)
silgk3 <- SIL.F(dataX.norm4z,gkpfcm.result3$u)
silgk4 <- SIL.F(dataX.norm4z,gkpfcm.result4$u)
silgk5 <- SIL.F(dataX.norm4z,gkpfcm.result5$u)
silgk6 <- SIL.F(dataX.norm4z,gkpfcm.result6$u)
silgk7 <- SIL.F(dataX.norm4z,gkpfcm.result7$u)
silgk8 <- SIL.F(dataX.norm4z,gkpfcm.result8$u)
silgk9 <- SIL.F(dataX.norm4z,gkpfcm.result9$u)
silgk10 <- SIL.F(dataX.norm4z,gkpfcm.result10$u)
silhouetteGK=cbind.data.frame(silgk2,silgk3,silgk4,silgk5,silgk6,silgk7,silgk8,silgk9,silgk10)
GKsil=c(0,silgk2,silgk3,silgk4,silgk5,silgk6,silgk7,silgk8,silgk9,silgk10)
GKsilhouette=data.frame(k=c(1,2,3,4,5,6,7,8,9,10),GKsil)
library(ggpubr)
library(ggplot2)
library(magrittr)
plotsilgk=ggline(GKsilhouette, x = "k", y = "GKsil",group = 1,color = "blue",
                 title = "Optimal number of cluster", xlab="Number of clusters k", ylab="Average silhouette width")
p <- plotsilgk + geom_vline(xintercept = 2, linetype=2, color = "blue")

###PseudoFGK#
datainputpsgk2=cbind(dataX.norm4z,gkpfcm.result2$cluster)
datainputpsgk3=cbind(dataX.norm4z,gkpfcm.result3$cluster)
datainputpsgk4=cbind(dataX.norm4z,gkpfcm.result4$cluster)
datainputpsgk5=cbind(dataX.norm4z,gkpfcm.result5$cluster)
datainputpsgk6=cbind(dataX.norm4z,gkpfcm.result6$cluster)
ujips.gk2 = icdrate(datainputpsgk2,154,2)
ujips.gk3 = icdrate(datainputpsgk3,154,3)
ujips.gk4 = icdrate(datainputpsgk4,154,4)
ujips.gk5 = icdrate(datainputpsgk5,154,5)
ujips.gk6 = icdrate(datainputpsgk6,154,6)
kOptim.gk=cbind.data.frame(kgk2=c(ujips.gk2$pseudof,ujips.gk2$icdrate),kgk3=c(ujips.gk3$pseudof,ujips.gk3$icdrate),kgk4=c(ujips.gk4$pseudof,ujips.gk4$icdrate),kgk5=c(ujips.gk5$pseudof,ujips.gk5$icdrate),kgk6=c(ujips.gk6$pseudof,ujips.gk6$icdrate))
rownames(kOptim.gk)=c("pseudof","icdrate")

###Plot Gustafson Kessel###
##Visualizing##
library(rgl)
library(car)
library(carData)
cnames=c("dim1","dim2","dim3","dim4")
ind = sort(apply(dataX.norm4z,2,var),index.return=TRUE, decreasing = TRUE)$ix[1:3]
labs = cnames[ind]
labels = gkpfcm.result6$cluster
sset = dataX.norm4z[,ind]
x = sset[,1]
y = sset[,2]
z = sset[,3]
par3d("windowRect"= c(0,0,400,400))
scatter3d(x = x, y = y, z = z, 
          xlab=labs[1], ylab=labs[2], zlab = labs[3],
          labels = NULL,
          groups = as.factor(labels), 
          surface = FALSE, 
          grid = FALSE,
          ellipsoid=TRUE)
rgl.snapshot(filename = "plot_GustafsonKessel6.png")

