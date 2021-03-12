# Import data and packages
library(raster)
library(ggplot2)
library(ggfortify)
library(psych)
library(xtable)
library(dplyr)
library(boot)
library(corrplot)
library(rgl)
library(robustbase)
library(rrcov)
library(cowplot)
library(parallel)
library(reshape2)

theme_set(theme_gray())

setwd("~/Multivariate/Project")

source("HelperFunc.R")

osa_alueetPoly <- shapefile(paste0("C:\\Users\\Mikko Kaivola\\Documents\\ArcGIS\\mapdata\\osa_alueet_tiedot",
                            "\\osa_alueet_merip\\JoinedMaps\\RemovePoliceUnknown\\OnlyWhatIsKnown\\FinalJoined2.shp"))

osaAlueData <- osa_alueetPoly@data

# Clean data

osaAlueDataClean <- osaAlueData %>%
  select(-c(1:6,8)) %>%
  select(-AECount, -TBSTCount, -Aream2, -Traf2017, -HHI2017)

osaAlueDataClean[,c(2:ncol(osaAlueDataClean))] <- lapply(osaAlueDataClean[,c(2:ncol(osaAlueDataClean))],
                                                         function(x) as.numeric(as.character(x)))

osaAlueTraining <- osaAlueDataClean %>%
  select(-contains("2016")) %>%
  select(-contains("2014")) %>%
  filter(GINI2015 != 0 )

#write.csv(dataMatrix, file = "HelsinkiRegionData.csv", row.names = F)

# Univariate and bivariate analysis ####

univData <- select (describe(osaAlueTraining[,-1]), -c(vars,n,trimmed,se,range))

tableFactory(univData,"univtable.tex")

# Bootstrap confidence intervals for univariate

cpuNum <- detectCores()
bootClust <- makeCluster(cpuNum)

bootstrapInts <- function(data,cpuNum) {
 bootObj <- boot(data, function(x,n) as.matrix( select (describe(x[n]), -c(vars,n,trimmed,se,range)
  ) )[1,],10000 , parallel = "snow", ncpus = cpuNum, cl = bootClust)
  sapply(1:8,function(x) boot.ci(bootObj,index = x, conf = 0.95, type = "perc")$percent[-c(1:3)])
}

clusterExport(bootClust,list("select","describe"))

confInts <- (apply(osaAlueTraining[,-1],2, 
      function(x) bootstrapInts(x, 4)))

stopCluster(bootClust)

confMatr <- t(apply(confInts,2,function(x) {
  paste("(",round(x[c(T,F)],2)," , ",round(x[c(F,T)],2),")", sep = "")
}))

colnames(confMatr) <- colnames(univData)

tableFactory(confMatr,"BootstrapConfUni.tex")

dotGene <- function(varName, label) {
  ggplot(data = osaAlueTraining, aes_string(x = "NIMI_ISO", y = varName)) + 
    geom_point(col = "black", size = 2) + 
    geom_segment(aes_string(x="NIMI_ISO", 
                   xend="NIMI_ISO", 
                   y=min(varName), 
                   yend=max(varName)), 
                linetype="dashed", 
               size=0.1) + ylab(label) + xlab("Area") +
    coord_flip()
}

AEdot <- dotGene("AEDensity","Alcohol establishment density")
AEhist <- fancyHist(osaAlueTraining$AEDensity,"Alcohol establishment density",10000/2)
plot_grid(AEdot,AEhist)
TBSTdot <- dotGene("TBSTDen","Station density")
TBSThist <- fancyHist(osaAlueTraining$TBSTDen,"Station density",10000*0.8)
plot_grid(TBSTdot,TBSThist)
Omaisdot <- dotGene("OmaisR2017","Property crimes")
Omaishist <- fancyHist(osaAlueTraining$OmaisR2017,"Number of property crimes",150)
plot_grid(Omaisdot,Omaishist)
Henkidot <- dotGene("HenkiR2017","Violent crimes")
Henkihist <- fancyHist(osaAlueTraining$HenkiR2017,"Violent crimes",20)
plot_grid(Henkidot,Henkihist)

# Bootstrap conf. int. for correlations (explanatory variables)

bootClust <- makeCluster(detectCores())

corConf <- apply(select(osaAlueTraining[,-1],-HenkiR2017,-OmaisR2017),2,function(x) {
  apply(select(osaAlueTraining[,-1],-HenkiR2017,-OmaisR2017),2, function(z) { res <- boot.ci(boot(cbind(x,z), 
                                                         function(y,n) cor(y[n,1],y[n,2], method = "pearson"), 10000,
                                                         parallel = "snow",
                                                         ncpus = 4,
                                                         cl = bootClust), conf = 0.95, type = "perc")$percent[-c(1:3)]
  if (is.null(res)) {return(c(1,1))} else {return(res)}
    })
})

stopCluster(bootClust)

corConfMat <- apply(corConf,2,function(x) {
  paste("(",round(x[c(T,F)],2)," , ",round(x[c(F,T)],2),")", sep = "")
})

xind <- row(corConfMat)
yind <- (ncol(corConfMat) + 1) - col(corConfMat)

corMat <- cor(select(osaAlueTraining[,-1],-HenkiR2017,-OmaisR2017), method = "spearman")

corrplot(corMat,method = "number")
text(xind,yind,corConfMat, pos = 1, cex = 0.6)

# Scatter plots

anomal1 <- ggplot(data=osaAlueTraining,aes(x = TBSTDen, y = Work2015)) + geom_point(size = 1.5, color = "blue", alpha = 0.6) + 
  xlab("Station density") + ylab("Number of jobs")

anomal2 <- ggplot(data=osaAlueTraining,aes(x = TBSTDen, y = AEDensity)) + geom_point(size = 1.5, color = "blue", alpha = 0.6) + 
  xlab("Station density") + ylab("Alcohol establishment density")

plot_grid(anomal1,anomal2)

# Conf. ints for correlations (response varib.)

bootClust <- makeCluster(cpuNum)

corResponse <- apply(select(osaAlueTraining,HenkiR2017,OmaisR2017),2,function(x) 
  apply(select(osaAlueTraining[,-1],-HenkiR2017,-OmaisR2017),2, function(z) boot.ci(boot(cbind(x,z), 
                                                                                                  function(y,n) cor(y[n,1],y[n,2], method = "pearson"), 10000,
                                                                                                  parallel = "snow",
                                                                                                  ncpus = 4,
                                                                                                  cl = bootClust), conf = 0.95, type = "perc")$percent[-c(1:3)]
  ))


stopCluster(bootClust)

corResponseMat <- apply(corResponse,2,function(x) {
  paste("(",round(x[c(T,F)],2)," , ",round(x[c(F,T)],2),")", sep = "")
})

rownames(corResponseMat) <- colnames(select(osaAlueTraining[,-1],-HenkiR2017,-OmaisR2017))

corResponseMat <- t(corResponseMat)

corRespMat <- t(apply(select(osaAlueTraining,HenkiR2017,OmaisR2017),2,function(x) 
  apply(select(osaAlueTraining[,-1],-HenkiR2017,-OmaisR2017),2, function(z) cor(x,z, method = "pearson"))))

yindR <- (nrow(corResponseMat) + 1) - row(corResponseMat)
xindR <- col(corResponseMat)

pdf("corResponse.pdf", height=3, width=10)
corrplot(corRespMat,method = "number")
text(xindR,yindR,corResponseMat, pos = 1, cex = 0.6)
dev.off()

YHhenki <- ggplot(osaAlueTraining,aes(x = YH2017, y = HenkiR2017)) + geom_point(size = 1.5, color = "blue", alpha = 0.6) + 
  xlab("Number of single parents") + ylab("Violent crimes")

AEHenki <- ggplot(osaAlueTraining,aes(x = AEDensity, y = HenkiR2017)) + geom_point(size = 1.5, color = "blue", alpha = 0.6) + 
  xlab("Alcohol establishment density") + ylab("Violent crimes")

plot_grid(YHhenki,AEHenki)

# ----

# Multivariate

trainingExplan <- select (osaAlueTraining[,-1], -HenkiR2017, -OmaisR2017)
rownames(trainingExplan) <- osaAlueTraining$NIMI_ISO

# PCA
?prcomp
pcaResult <- prcomp(trainingExplan,scale. = T)
summary(pcaResult)
varianceTable <- rbind(pcaResult$sdev,pcaResult$sdev^2/sum(pcaResult$sdev^2),
                       cumsum(pcaResult$sdev^2/sum(pcaResult$sdev^2)))
colnames(varianceTable) <- paste("PC",1:10, sep = "")
rownames(varianceTable) <- c("Standard deviation", "Proportion of variance",
                             "Cumulative proportion")
tableFactory(varianceTable,"PCAvarTable.tex")

loadingS <- sweep(pcaResult$rotation,2,pcaResult$sdev,FUN = "*")
corrplot(loadingS,method = "number")

varianceAcc <- t(apply(loadingS^2,1,cumsum))
corrplot(varianceAcc,is.corr = F,method = "number")
biplot(pcaResult)

contribInd <- sweep((pcaResult$x)^2 * 1/(nrow(osaAlueTraining) - 1),2,pcaResult$sdev^2,FUN = "/")

contrFunc <- function(index) {
sortedContr <- sort(contribInd[,index], decreasing = T)[1:10]

ggplot(data.frame(contri = sortedContr, names = names(sortedContr)), aes(x = reorder(names,-contri), y = contri)) + 
  geom_bar(stat = "identity", fill = "blue", alpha = 0.8) + xlab("Observational units") + ylab("Contribution") +
  geom_hline(yintercept = 1/nrow(contribInd), color = "red") + labs(title = paste0("PC",index)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylim(0,1)
}

plot_grid(contrFunc(1),contrFunc(2),contrFunc(3),contrFunc(4))

qualityInd <- sweep((pcaResult$x)^2, 1, rowSums((pcaResult$x)^2), FUN = "/")
firstCompQuality <- data.frame(qualityInd[,1:4], names = rownames(qualityInd))
qualityMelted <- melt(firstCompQuality,id.vars = "names")
colnames(qualityMelted)[2] <- "Component"

ggplot(qualityMelted, aes(x = reorder(names,-value), y = value, fill = Component)) + geom_bar(stat = "identity",alpha = 0.8) + 
  xlab("Observational units") + ylab("Quality of representation") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) # + 
  #scale_fill_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3"))

#ggplot(firstCompQuality, aes(x = reorder(names,-Quality), y = Quality)) + 
#  geom_bar(stat = "identity", fill = "blue", alpha = 0.8) + 
#  xlab("Observational units") + ylab("Quality of representation") + 
#  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(data = data.frame(pcaResult$x), aes(x = PC1, y = PC4, label = rownames(trainingExplan))) + geom_point() +
  geom_text() + xlim(-8.2,4.5)

open3d()
plot3d(pcaResult$x[,1:3])
identify3d(pcaResult$scores[,1:3], labels = osaAlueTraining$NIMI_ISO)
 
autoplot(pcaResult,loadings = T, loadings.colour = 'blue', loadings.label = T, 
         loadings.label.repel = T, x = 1, y = 4)

# Confints for response correlation with scores

bootClust <- makeCluster(detectCores())


scores <- pcaResult$x

corResponse <- apply(select(osaAlueTraining,HenkiR2017,OmaisR2017),2,function(x) 
  apply(scores,2, function(z) boot.ci(boot(cbind(x,z), 
                                                                                         function(y,n) cor(y[n,1],y[n,2], method = "pearson"), 10000,
                                                                                         parallel = "snow",
                                                                                         ncpus = 4,
                                                                                         cl = bootClust), conf = 0.95, type = "perc")$percent[-c(1:3)]
  ))


stopCluster(bootClust)

corResponseMat <- apply(corResponse,2,function(x) {
  paste("(",round(x[c(T,F)],2)," , ",round(x[c(F,T)],2),")", sep = "")
})

rownames(corResponseMat) <- colnames(scores)

corResponseMat <- t(corResponseMat)

corRespMat <- t(apply(select(osaAlueTraining,HenkiR2017,OmaisR2017),2,function(x) 
  apply(scores,2, function(z) cor(x,z, method = "pearson"))))

yindR <- (nrow(corResponseMat) + 1) - row(corResponseMat)
xindR <- col(corResponseMat)

pdf("corResponsePCA.pdf", height=3, width=10)
corrplot(corRespMat,method = "number")
text(xindR,yindR,corResponseMat, pos = 1, cex = 0.6)
dev.off()

p4 <- ggplot(data.frame(scores,Omais2017 = osaAlueTraining$OmaisR2017), aes(x = PC4, y = Omais2017 )) + 
  geom_point() + ylab("Property crimes")

plot_grid(p1,p2,p3,p4)

p4
# Robust tests ####
pcaRobust <- PcaCov(trainingExplan ,scale=T,cov.control = CovControlMcd(alpha=1/2))
summary(pcaRobust)
biplot(pcaRobust)
loadingSR <- sweep(pcaRobust$loadings,2,sqrt(pcaRobust$eigenvalues), FUN ="*")
corrplot(loadingSR,method = "number")
# ----

# Regression ####

scoreValues <- data.frame(pcaResult$x)
normalizedCrime <- data.frame((select (osaAlueTraining[,-1], HenkiR2017, OmaisR2017)))
rownames(normalizedCrime) <- osaAlueTraining$NIMI_ISO

regresData <- cbind(scoreValues,normalizedCrime)
corrplot(cor(regresData,method = "pearson"),method = "number")
corrplot(covMcd(regresData,cor = T)$cor, method = "number")

model <- ltsReg(HenkiR2017~ PC1 + PC2 + PC3 + PC4 , data = regresData)
summary(model)
plot(fitted(model),resid(model))
# ----
