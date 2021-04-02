############################################################################
###set the path
getwd()
setwd("E:/targetprofiling")
rm(list = ls(all = TRUE))


path = "E:/data/20170612chembl_bindingdb_all/CHEMBL_ BINDINDDB/chembl_bindingdb_maccs/bindingdb_chembl_unis_all"


pid = read.csv("E:/targetprofiling/proteinid.csv", sep = ",",header = TRUE)
pid = as.vector(pid$uniprotID)


require(gbm)
############################################################################
GBMmaccsper = list()

for (i in pid){
  dat = read.csv(paste0(path,"/",i,".txt"), sep = "\t",header = TRUE,stringsAsFactors=FALSE)
  Xdata = dat[,-c(1:7)]
  temp = as.numeric(dat[,4])
  ydata = (temp < median(temp))
  X = data.frame(y = ydata,Xdata)
  GBM = gbm(y~.,data = X,distribution = "bernoulli",cv.folds = 5,n.trees = 5000,shrinkage=0.001,bag.fraction =0.5,n.cores=7)
  save(GBM,file =paste0("E:/targetprofiling/model/GBM","/",i,".maccs.Rdata"))
  ypre = sign(GBM$cv.fitted)
  ypre[ypre==-1] = 0
  GBMmaccsper[i] = sum((ypre-ydata)!=0)/length(ydata)
  print(i)
}

############################################################################

save(GBMmaccsper,file ="E:/targetprofiling/model/GBMmaccsper.Rdata")
