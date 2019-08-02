#!/user/bin/Rscript
library(tidyverse)
library(minpack.lm)
library(dplyr)
library(jsonlite)
library(arrayhelpers)
library(abind)
#library(rlang)
##counts<- read.csv("./Example.csv", fileEncoding = "UTF-8-BOM")
##Load and derive information from JSON
args<- commandArgs(TRUE)
json<-args[1]
groupID<-args[2]
#jsondf<-fromJSON("./sample.json")
jsondf<-fromJSON(json)
compounds<-jsondf$compounds
compoundorder<-cbind(compounds, c(1:length(compounds)))
data_replicates<-jsondf$replicates
data_outliers<-jsondf$outliers
col_num<-ncol(data_replicates)
row_num<-nrow(data_replicates)
rep_num<-dim(data_replicates)[3]
#Manualplate indicator, 0 by default, gets set to 1 if the data is a manual plate.
Manualplate<-0

Manualdata<-function (manualjsondf){
  data_replicates_manual<-array(dim = c(3,12,3))
  for (i in 1:length(manualjsondf)){
    for (j in 1:nrow((manualjsondf[[i]]))){
      for(k in 1:ncol((manualjsondf[[i]]))){
        data_replicates_manual[k,j,i]<-(as.numeric(manualjsondf[[i]][j,k]))
      }
    }
  }
  return(data_replicates_manual)
}

if(is.null(rep_num)==TRUE){
  data_replicates<-Manualdata(data_replicates)
  data_outliers<-Manualdata(data_outliers)
  rep_num<-dim(data_replicates)[3]
  col_num<-ncol(data_replicates)
  row_num<-nrow(data_replicates)
  Manualplate<-1
}

##Create empty array to check for outliers with. Array is either 3D or 2D, but will match data from JSON
data_replicates_checked<-array(NA, dim=dim(data_replicates))
##Check for outliers on the outlier grid. If it is an outlier, the data grid has the corresponding value
##set to NA
outliercheck<- function(replicates, outliers){
  #if (rep_num>0){
  outliers<-na.omit(outliers)
  for (i in 1:dim(replicates)[1]) {
    for(j in 1:dim(replicates)[2]) {
      for (k in 1:dim(replicates)[3]){
      if(is.na(outliers[i,j,k])==TRUE){
        data_replicates_checked[i,j,k]<-NA
      }
      else if(outliers[i,j,k] == 0){
        data_replicates_checked[i,j,k]<-as.numeric(as.character(replicates[i,j,k]))
      }
      else if(outliers[i,j,k] !=0){
        data_replicates_checked[i,j,k]<-NA
      }
      }
    }
  }
  return(data_replicates_checked)
  }
  # else if (is.na(rep_num) == TRUE){
  #   for (i in 1:row_num) {
  #     for(j in 1:col_num) {
  #         if(outliers[i,j] == 0){
  #           data_replicates_checked[i,j]<-as.numeric(as.character(replicates[i,j]))}
  #         else{
  #           data_replicates_checked[i,j]<-NA
  #         }
  #     }
  #     return(data_replicates_checked)
  # }
  # }
data_replicates_checked<-outliercheck(data_replicates, data_outliers)
#With outliers set as NA, bind the compound number to the first column in each plate and replace with numbered
#order 1-8. This differs for manual plates, so we check the indicator and decide.
if(Manualplate == 0){
  data_replicates_cmpd<-abind(array(compoundorder[,2], replace(dim(data_replicates_checked),2,1)), 
  data_replicates_checked, along = 2)
}else if(Manualplate == 1){
  data_replicates_cmpd<-abind(array(rep(1:nrow(compoundorder), each = 3), dim=c(3,1,rep_num)),
  data_replicates_checked, along = 2)
  }
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
#Convert the data columns into numeric and create average values
data_replicates[,2:ncol(data_replicates)]<-lapply(data_replicates[,2:ncol(data_replicates)], function(x) as.numeric(as.character(x)))
data_avgcounts<-data_replicates %>% group_by(V1)%>%summarise_all(mean, na.rm=TRUE)%>%
  arrange(match(V1,c(compoundorder[,2])), desc(V1))
#Prep assay information from JSON, most importantly the constant used in the equation
curies<-as.numeric(jsondf$hot_count)/(2.22*10^12)
mmols<-curies/as.numeric(jsondf$hot_activity)
hotnM<-mmols*10^6/(as.numeric(jsondf$hot_volume)*10^-6)
#Constants to be used in logKi equation
Kd<-as.numeric(jsondf$dissociation_constant)
kvalue<-log10(1+(hotnM/Kd))
#Prep concentration/X-values from JSON
concentrations<-as.data.frame(array(as.numeric(jsondf$concentrations), dim = dim(jsondf$concentrations)))
concentrations<-cbind(compoundorder[,2], concentrations)
#Check and set concentrations to NA if the corresponding Y value is NA
for(i in 1:nrow(data_avgcounts)){
  for(j in 1:ncol(data_avgcounts)){
    if(is.na(data_avgcounts[i,j]== TRUE)){
      concentrations[i,j]<-NaN
    }
  }
}
#Global Regression requires stacking Y values and X values
#Vectorize all Y values, drop NaN or NA values, and match them with the same X-values. Drop X-values that
#correspond to a dropped Y-value.
Yvals<-setNames(split(data_avgcounts, f = data_avgcounts$V1),
                  paste0("Y",1:length(compounds)))
Xvals<-setNames(split(concentrations, f = data_avgcounts$V1),
                        paste0("X",1:length(compounds)))
for(i in 1:length(compounds)){
  Yvals[[i]]<-Filter(function(x)!all(is.na(x)), Yvals[[i]])
  Xvals[[i]]<-Filter(function(x)!all(is.na(x)), Xvals[[i]])
}
#Stack all average counts into a single Y vector, same with concentrations into a single X vector
Yavg<-c()
Xavg<-c()
for(i in 1:length(compounds)){
  insertY<-unlist(Yvals[[i]][2:length(Yvals[[i]])])
  Yavg<-c(Yavg,insertY)
  insertX<-unlist(Xvals[[i]][2:length(Xvals[[i]])])
  Xavg<-c(Xavg,insertX)
}
#Yavg<c(t(data_replicates)
#Xavg<-c(t(concentrations))
#Global regression will require indicator matrices to pair individual Ki's to 
#their respective data sets. However, for shared parameters this is not needed.
#We make a giant list where each vector is going to be a single indicator vector of 1's and 0's
#These vectors are similar to design matrices found in linear regression/ANOVA
#1 indicates this data set value is for this parameter, while 0 indicates this value is not for this parameter.
Indicatorlist<-list()
for (i in 1:length(compounds)){
  Indicatorlist[[i]]<-rep(0, length(Yavg))
}
index<-0
for (i in 1:length(compounds)){
  indexprev<-index+1
  index<-index+(length(Yvals[[i]][2,])-1)
  Indicatorlist[[i]][indexprev:index]<-1
}
#Indicatortable<-matrix(0, nrow = length(compounds), ncol = 12*length(compounds))
#for(i in 1:length(compounds)){
#  Indicatortable[i, (12*i-11):(12*i)]<-1
#}
#Generate Logki Parameter table. First set all ki's to -7. Then add top and bottom estimates
kiparams<-c()
startingki<-c()
for(i in 1:length(compounds)){
  kiparams[i]<-cbind(paste0("logKi", i))
  startingki[i]<-cbind(-7)
}
paramslist<-as.list(setNames(startingki, kiparams))
bot_est<-0
top_est<-mean(data_avgcounts$V2, na.rm = TRUE)
paramslist[["Top"]]<-top_est
paramslist[["Bottom"]]<-bot_est
#Now we need to insert the correct number of variables into our formula. A huge pain.
# paramsindics<-c()
# for (i in 1:length(compounds)){
#   paramsindics[i]<-c(paste0("params$logKi",i,"*","Indicatorlist[[",i,",]]"))
# }
#Now we can start our predictions. There can be anywhere form 2-7 compounds so we need to create an equation
#for each situation. There is probably a way to do this programmatically but I have decided to use
#a brute force method for now. Bascially, a large if-then loop to match the compounds to the equation.
if (length(compounds) == 8){
  getPred<- function(params, xx) {
    (params$Bottom) + ((params$Top-params$Bottom)/
                         (1+(10^(xx-params$logKi1*Indicatorlist[[1]]-params$logKi2*Indicatorlist[[2]]
                                 -params$logKi3*Indicatorlist[[3]]
                                 -params$logKi4*Indicatorlist[[4]]
                                 -params$logKi5*Indicatorlist[[5]]
                                 -params$logKi6*Indicatorlist[[6]]
                                 -params$logKi7*Indicatorlist[[7]]
                                 -params$logKi8*Indicatorlist[[8]]-kvalue))))
  }

  } else if (length(compounds) == 7){
    getPred<- function(params, xx) {
      (params$Bottom) + ((params$Top-params$Bottom)/
                           (1+(10^(xx-params$logKi1*Indicatorlist[[1]]-params$logKi2*Indicatorlist[[2]]
                                   -params$logKi3*Indicatorlist[[3]]
                                   -params$logKi4*Indicatorlist[[4]]
                                   -params$logKi5*Indicatorlist[[5]]
                                   -params$logKi6*Indicatorlist[[6]]
                                   -params$logKi7*Indicatorlist[[7]]
                                   -kvalue))))
    }
  
   } else if (length(compounds) == 6){
      getPred<- function(params, xx) {
        (params$Bottom) + ((params$Top-params$Bottom)/
                             (1+(10^(xx-params$logKi1*Indicatorlist[[1]]-params$logKi2*Indicatorlist[[2]]
                                     -params$logKi3*Indicatorlist[[3]]
                                     -params$logKi4*Indicatorlist[[4]]
                                     -params$logKi5*Indicatorlist[[5]]
                                     -params$logKi6*Indicatorlist[[6]]
                                     -kvalue))))
      }
    
  } else if (length(compounds) == 5){
    getPred<- function(params, xx) {
      (params$Bottom) + ((params$Top-params$Bottom)/
                           (1+(10^(xx-params$logKi1*Indicatorlist[[1]]-params$logKi2*Indicatorlist[[2]]
                                   -params$logKi3*Indicatorlist[[3]]
                                   -params$logKi4*Indicatorlist[[4]]
                                   -params$logKi5*Indicatorlist[[5]]
                                   -kvalue))))
    }
  
  } else if (length(compounds) == 4){
    getPred<- function(params, xx) {
      (params$Bottom) + ((params$Top-params$Bottom)/
                           (1+(10^(xx-params$logKi1*Indicatorlist[[1]]-params$logKi2*Indicatorlist[[2]]
                                   -params$logKi3*Indicatorlist[[3]]
                                   -params$logKi4*Indicatorlist[[4]]
                                   -kvalue))))
    }
  
  } else if (length(compounds) == 3){
    getPred<- function(params, xx) {
      (params$Bottom) + ((params$Top-params$Bottom)/
                           (1+(10^(xx-params$logKi1*Indicatorlist[[1]]-params$logKi2*Indicatorlist[[2]]
                                   -params$logKi3*Indicatorlist[[3]]
                                   -kvalue))))
    }
  
  } else if (length(compounds) == 2){
    getPred<- function(params, xx) {
      (params$Bottom) + ((params$Top-params$Bottom)/
                           (1+(10^(xx-params$logKi1*Indicatorlist[[1]]-params$logKi2*Indicatorlist[[2]]
                                   -kvalue))))
    }
  
  }else {
    print("ERROR IN EQUATION GENERATION")
  }


# subtractionf<-(eval(as.character(paste(paramsindics, collapse="-"))))
# as.expression((quote(params$Bottom)) + ((quote(params$Top-params$Bottom))/
#                        (1+(10^(quote(xx)-subtractionf-quote(kvalue))))))

#Create simulated values based on equation and starting estimates.
simDnoisy <- getPred(paramslist, Xavg)
#The residual function used to gauage our residuals
residFun<- function(p, observed, xx) {
  observed - getPred(p,xx)
}
#The actual regression function. From the nls.lm package
nls.out<- nls.lm(paramslist, fn = residFun, observed = Yavg, xx = Xavg)
#summary(nls.out)
#Grab predicted Y-values using the updated parameters of the equation.
getPredsingle<- function(params, xx) {
  (params$Bottom) + ((params$Top-params$Bottom)/
  (1+(10^(xx-params$logKi-kvalue))))
}
#Compile results.
paramslisttotal<-list()
for(i in 1:length(compounds)){
  paramslisttotal[[i]]<-list(logKi=coef(nls.out)[[paste0("logKi",i)]],Top=coef(nls.out)[["Top"]],Bottom=coef(nls.out)[["Bottom"]])
}
output<-list()
spl_fxns<-list()
output_assembly<-paramslisttotal
#Rename output_assembly using compound names. Now we can construct our JSON to be passed back.
names(output_assembly)<-c(compounds)
for(i in 1:length(compounds)){
output[[i]]<-getPredsingle(paramslisttotal[[i]], unlist(Xvals[[i]][,2:length(Xvals[[i]])]))
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
  finaloutput[[6]][[i]]<-as.vector(as.numeric(Xvals[[i]][,2:length(Xvals[[i]])]))
}
#Rename our output
names(finaloutput)<-c("Compounds", "logKi","Top","Bottom","YValues","XValues")
jsonoutput<-toJSON(finaloutput, pretty=TRUE, auto_unbox = TRUE)
write(jsonoutput, paste0(paste0("/path/to/save/",groupID),".json"))
#p<-ggplot(data.frame(x=x, y=Yavg), mapping = aes(Xavg, Yavg)) + 
  #geom_point()
#for (i in 1:length(compounds)){
  #p<-p+stat_function(fun = spl_fxns[[i]])
#}

