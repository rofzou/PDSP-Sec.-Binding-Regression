#!/user/bin/Rscript
library(tidyverse)
library(minpack.lm)
library(dplyr)
library(jsonlite)
library(arrayhelpers)
library(abind)
##counts<- read.csv("./Example.csv", fileEncoding = "UTF-8-BOM")
##Load and derive information from JSON
args<- commandArgs(TRUE)
json<-args[1]
groupID<-args[2]
#jsondf<-fromJSON("./sample.json")
jsondf<-fromJSON(json)
compounds<-jsondf$compounds
data_replicates<-jsondf$replicates
data_outliers<-jsondf$outliers
col_num<-ncol(data_replicates)
row_num<-nrow(data_replicates)
rep_num<-dim(data_replicates)[3]
##Create empty array to check for outliers with. Array is either 3D or 2D, but will match data from JSON
data_replicates_checked<-array(NA, dim=dim(data_replicates))
outliercheck<- function(replicates, outliers){
  if (rep_num>0){
  for (i in 1:row_num) {
    for(j in 1:col_num) {
      for (k in 1:rep_num){
      if(outliers[i,j,k] == 0){
        data_replicates_checked[i,j,k]<-as.numeric(as.character(replicates[i,j,k]))}
      else{
        data_replicates_checked[i,j,k]<-NA
      }
    }
    }
  }
    return(data_replicates_checked)
  }
  else if (is.na(rep_num) == TRUE){
    for (i in 1:row_num) {
      for(j in 1:col_num) {
          if(outliers[i,j] == 0){
            data_replicates_checked[i,j]<-as.numeric(as.character(replicates[i,j]))}
          else{
            data_replicates_checked[i,j]<-NA
          }
      }
      return(data_replicates_checked)
  }
  }
}
data_replicates_checked<-outliercheck(data_replicates, data_outliers)
#With outliers set as NA, bind the compound number to the first column in each plate
data_replicates_cmpd<-abind(array(compounds, replace(dim(data_replicates_checked),2,1)), 
                       data_replicates_checked, along = 2)
#Arrange the data into a single 2D array depeding on 3D-ness of original array
data_replicates<-data_replicates_cmpd[,,1]
if(rep_num>0){
for (k in 2:rep_num){
  data_replicates<-rbind(data_replicates, data_replicates_cmpd[,,k])
}
}
if (is.na(rep_num)== TRUE){
  data_replicates<-data_replicates_cmpd[,,1]
}
data_replicates<-as.data.frame(data_replicates)
#Convert the data columns into numeric and create average values and standard deviation
data_replicates[,2:13]<-lapply(data_replicates[,2:13], function(x) as.numeric(as.character(x)))
data_avgcounts<-data_replicates %>% group_by(V1) %>% summarise_all(mean, na.rm=TRUE)
#Prep assay information from JSON
curies<-as.numeric(jsondf$hot_count)/(2.22*10^12)
mmols<-curies/as.numeric(jsondf$hot_activity)
hotnM<-mmols*10^6/(as.numeric(jsondf$hot_volume)*10^-6)
#Constant to be used in logKi equation
Kd<-as.numeric(jsondf$dissociation_constant)
kvalue<-log(1+(hotnM/Kd))
#Prep concentration/X-values from JSON
concentrations<-as.data.frame(array(as.numeric(jsondf$concentrations), dim = dim(jsondf$concentrations)))
#Global Regression requires stacking Y values and X values
#Stack all average counts into a single Y column, same with concentrations into a single X column
Yavg<-c(t(data_avgcounts[2:13]))
Xavg<-c(t(concentrations))
#Global regression will require indicator matrices to pair individual Ki's to 
#their respective data sets. However, for shared parameters this is not needed.
#We make a giant matrix where each row is going to be a single indicator matrix
#These matrices are similar to design matrices found in linear regression/ANOVA
Indicatortable<-matrix(0, nrow = length(compounds), ncol = 12*length(compounds))
for(i in 1:length(compounds)){
  Indicatortable[i, (12*i-11):(12*i)]<-1
}
#Generate Logki Parameter table. First set all ki's to -7. Then add top and bottom
kiparams<-c()
startingki<-c()
for(i in 1:length(compounds)){
  kiparams[i]<-cbind(paste0("logKi", i))
  startingki[i]<-cbind(-7)
}
paramslist<-as.list(setNames(startingki, kiparams))
bot_est<-0
top_est<-mean(data_avgcounts$V2)
paramslist[["Top"]]<-top_est
paramslist[["Bottom"]]<-bot_est
#Now we need to insert the correct number of variables into our formula. A huge pain.
paramsindics<-c()
for (i in 1:length(compounds)){
  paramsindics[i]<-c(paste0("params$logKi",i,"*","Indicatortable[",i,",]"))
}
#Now we can start our predictions. Write the first equation.All 8 parameters.
getPred<- function(params, xx) {
  (params$Bottom) + ((params$Top-params$Bottom)/
  (1+(10^(xx-params$logKi1*Indicatortable[1,]-params$logKi2*Indicatortable[2,]-params$logKi3*Indicatortable[3,]
                                                        -params$logKi4*Indicatortable[4,]-params$logKi5*Indicatortable[5,]
                                                        -params$logKi6*Indicatortable[6,]-
                                                        params$logKi7*Indicatortable[7,]
                                                        -params$logKi8*Indicatortable[8,]-kvalue))))
}
#getPred2<- function(params, xx) {
#  (params$Bottom) + ((params$Top-params$Bottom)/
                       #(1+(10^(xx-as.character(paste(paramsindics, collapse="-"))-kvalue))))
#}

simDnoisy <- getPred(paramslist, Xavg)
residFun<- function(p, observed, xx) {
  observed - getPred(p,xx)
}
nls.out<- nls.lm(paramslist, fn = residFun, observed = Yavg, xx = Xavg)
#summary(nls.out)
#Separate all data
getPredsingle<- function(params, xx) {
  (params$Bottom) + ((params$Top-params$Bottom)/
  (1+(10^(xx-params$logKi-kvalue))))
}
paramslisttotal<-list()
for(i in 1:length(compounds)){
  paramslisttotal[[i]]<-list(logKi=coef(nls.out)[[paste0("logKi",i)]],Top=coef(nls.out)[["Top"]],Bottom=coef(nls.out)[["Bottom"]])
}
output<-list()
spl_fxns<-list()
output_assembly<-paramslisttotal
names(output_assembly)<-c(compounds)
for(i in 1:length(compounds)){
output[[i]]<-getPredsingle(paramslisttotal[[i]], concentrations[i,])
spl_fxns[[i]] <- splinefun(concentrations[i,], output[[i]])
}
for(i in 1:length(compounds)){
  output_assembly[[i]][["Pred Y-Value"]]<-output[[i]]
}
finaloutput<-vector("list",6)
finaloutput[[1]]<-c(compounds)
finaloutput[[2]]<-as.vector(unlist(sapply(output_assembly,"[","logKi"),use.names = FALSE))
finaloutput[[3]]<-output_assembly[[1]][["Top"]]
finaloutput[[4]]<-output_assembly[[1]][["Bottom"]]
predictedyvalues<-sapply(output_assembly,"[","Pred Y-Value")
for(i in 1:length(compounds)){
  finaloutput[[5]][[i]]<-as.numeric(unlist(predictedyvalues[[i]]))
}
#finaloutput[[5]]<-as.numeric(sapply(output_assembly,"[","Pred Y-Value"))
for(i in 1:length(compounds)){
  finaloutput[[6]][[i]]<-as.vector(as.numeric(concentrations[i,]))
}
names(finaloutput)<-c("Compounds", "logKi","Top","Bottom","Y-Values","X-values")
#jsonoutput<-toJSON(output_assembly, pretty=TRUE, auto_unbox = TRUE)
jsonoutput<-toJSON(finaloutput, pretty=TRUE, auto_unbox = TRUE)
write(jsonoutput, paste0(paste0("/path/to/save/",groupID),".json"))
#p<-ggplot(data.frame(x=x, y=Yavg), mapping = aes(Xavg, Yavg)) + 
  #geom_point()
#for (i in 1:length(compounds)){
  #p<-p+stat_function(fun = spl_fxns[[i]])
#}

