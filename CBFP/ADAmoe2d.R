############################################################################
###set the path
getwd()
setwd("E:/targetprofiling")
rm(list = ls(all = TRUE))


path = "E:/data/20170612chembl_bindingdb_all/CHEMBL_ BINDINDDB/chembl_bindingdb_moe2d/bindingdb_chembl_unis_all"


pid = read.csv("E:/targetprofiling/proteinid.csv", sep = ",",header = TRUE)
pid = as.vector(pid$uniprotID)


require(gbm)
############################################################################
ADAmoe2dper = list()

for (i in pid){
  dat = read.csv(paste0(path,"/",i,".txt"), sep = "\t",header = TRUE,stringsAsFactors=FALSE)
  Xdata = dat[,-c(1:7)]
  temp = as.numeric(dat[,4])
  ydata = (temp < median(temp))
  rnanid = which(is.na(Xdata),arr.ind=TRUE)
  if (dim(rnanid)[1]>1){
    Xdata = Xdata[-rnanid[,1],]
    ydata = ydata[-rnanid[,1]]
  }
  X = data.frame(y = ydata,Xdata)
  GBM = gbm(y~.,data = X,distribution = "adaboost",cv.folds = 5,n.trees = 5000,shrinkage=0.001,bag.fraction =0.5,n.cores=7)
  save(GBM,file =paste0("E:/targetprofiling/model/ADA","/",i,".moe2d.Rdata"))
  ypre = sign(GBM$cv.fitted)
  ypre[ypre==-1] = 0
  ADAmoe2dper[i] = sum((ypre-ydata)!=0)/length(ydata)
  print(i)
}

############################################################################

save(ADAmoe2dper,file ="E:/targetprofiling/model/ADAmoe2dper.Rdata")