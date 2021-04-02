############################################################################
###set the path
getwd()
setwd("E:/targetprofiling/cbf")
rm(list = ls(all = TRUE))

#####################################################
##load the data 


data = read.csv("E:/targetprofiling/cbf/data/descriptors/specs_cansmiles_wash_moe2d.txt", sep = "\t",header = TRUE, stringsAsFactors = FALSE)
data1 = as.matrix(data[,-c(1:3)])

data = read.csv("E:/targetprofiling/cbf/data/descriptors/specs_cansmiles_wash_cats.txt", sep = "\t",header = TRUE, stringsAsFactors = FALSE)
data2 = as.matrix(data[,-c(1:3)])


data = read.csv("E:/targetprofiling/cbf/data/descriptors/specs_cansmiles_wash_maccs.txt", sep = "\t",header = TRUE, stringsAsFactors = FALSE)
data3 = as.matrix(data[,-c(1:3)])

drugsmi = data[,c(1:3)]

n_input = dim(drugsmi)[1]   ##the number of the samples predicted by our methods
####################################################
##load the protein ID
pid = read.csv("E:/targetprofiling/proteinid.csv", sep = ",",header = TRUE)
pid = as.vector(pid$uniprotID)


####################################################
library(randomForest)
####################################################
path = "E:/targetprofiling/model/RF/"
rf2d = matrix(0,n_input,832)

for (i in 1:length(pid)){
  load(paste0(path,pid[i],".moe2d.Rdata"))
  rf2d[,i] = predict(RF,data1,type = "prob")[,2]
  rm(RF)
}

####################################################
path = "E:/targetprofiling/model/RF/"
rfcats = matrix(0,n_input,832)

for (i in 1:length(pid)){
  load(paste0(path,pid[i],".cats.Rdata"))
  rfcats[,i] = predict(RF,data2,type = "prob")[,2]
  rm(RF)
}

####################################################
path = "E:/targetprofiling/model/RF/"
rfmaccs = matrix(0,n_input,832)

for (i in 1:length(pid)){
  load(paste0(path,pid[i],".maccs.Rdata"))
  rfmaccs[,i] = predict(RF,data3,type = "prob")[,2]
  rm(RF)
}
####################################################
library(gbm)
####################################################
path = "E:/targetprofiling/model/GBM/"
gbm2d = matrix(0,n_input,832)

for (i in 1:length(pid)){
  load(paste0(path,pid[i],".moe2d.Rdata"))
  best.iter <- gbm.perf(GBM,method="cv",plot.it = FALSE)
  gbm2d[,i] = predict(GBM,data.frame(data1),best.iter,type = "response")
  rm(GBM)
}

####################################################
path = "E:/targetprofiling/model/GBM/"
gbmcats = matrix(0,n_input,832)

for (i in 1:length(pid)){
  load(paste0(path,pid[i],".cats.Rdata"))
  best.iter <- gbm.perf(GBM,method="cv",plot.it = FALSE)
  gbmcats[,i] = predict(GBM,data.frame(data2),best.iter,type = "response")
  rm(GBM)
}

####################################################
path = "E:/targetprofiling/model/GBM/"
gbmmaccs = matrix(0,n_input,832)

for (i in 1:length(pid)){
  load(paste0(path,pid[i],".maccs.Rdata"))
  best.iter <- gbm.perf(GBM,method="cv",plot.it = FALSE)
  gbmmaccs[,i] = predict(GBM,data.frame(data3),best.iter,type = "response")
  rm(GBM)
}
####################################################

####################################################
path = "E:/targetprofiling/model/ADA/"
ada2d = matrix(0,n_input,832)

for (i in 1:length(pid)){
  load(paste0(path,pid[i],".moe2d.Rdata"))
  best.iter <- gbm.perf(GBM,method="cv",plot.it = FALSE)
  ada2d[,i] = predict(GBM,data.frame(data1),best.iter,type = "response")
  rm(GBM)
}

####################################################
path = "E:/targetprofiling/model/ADA/"
adacats = matrix(0,n_input,832)

for (i in 1:length(pid)){
  load(paste0(path,pid[i],".cats.Rdata"))
  best.iter <- gbm.perf(GBM,method="cv",plot.it = FALSE)
  adacats[,i] = predict(GBM,data.frame(data2),best.iter,type = "response")
  rm(GBM)
}

####################################################
path = "E:/targetprofiling/model/ADA/"
adamaccs = matrix(0,n_input,832)

for (i in 1:length(pid)){
  load(paste0(path,pid[i],".maccs.Rdata"))
  best.iter <- gbm.perf(GBM,method="cv",plot.it = FALSE)
  adamaccs[,i] = predict(GBM,data.frame(data3),best.iter,type = "response")
  rm(GBM)
}
####################################################

save(pid,drugsmi,rf2d,rfcats,rfmaccs,gbm2d,gbmcats,gbmmaccs,ada2d,adacats,adamaccs, file = "specs.RData")
