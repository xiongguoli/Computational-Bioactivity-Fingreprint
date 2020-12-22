############################################################################
###set the path
getwd()
setwd("E:/targetprofiling")
rm(list = ls(all = TRUE))


path = "E:/data/20170612chembl_bindingdb_all/CHEMBL_ BINDINDDB/chembl_bindingdb_cats/bindingdb_chembl_unis_all"


pid = read.csv("E:/targetprofiling/proteinid.csv", sep = ",",header = TRUE)
pid = as.vector(pid$uniprotID)


require(randomForest)
############################################################################
RFcatsper = list()

for (i in pid){
  dat = read.csv(paste0(path,"/",i,".txt"), sep = "\t",header = TRUE,stringsAsFactors=FALSE)
  Xdata = dat[,-c(1:7)]
  temp = as.numeric(dat[,4])
  ydata = as.factor(temp < median(temp))
  RF = randomForest(Xdata,ydata,ntree = 500, mtry = 30)
  save(RF,file =paste0("E:/targetprofiling/model/RF","/",i,".cats.Rdata"))
  RFcatsper[i] = RF$err.rate[500,1]
  print(i)
}

############################################################################

save(RFcatsper,file ="E:/targetprofiling/model/RFcatsper.Rdata")