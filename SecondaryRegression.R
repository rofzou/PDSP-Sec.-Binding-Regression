#!/user/bin/Rscript
library(tidyverse)
library(minpack.lm)
library(dplyr)
library(jsonlite)
library(arrayhelpers)
library(abind)
library(rlang)
##counts<- read.csv("./Example.csv", fileEncoding = "UTF-8-BOM")
##Load and derive information from JSON
args<- commandArgs(TRUE)
json<-args[1]
groupID<-args[2]
#This line loads the information manually so you can test it.
#jsondf<-fromJSON("./sample.json")
#We start extracting information from the json here.
jsondf<-fromJSON(json)
compounds<-jsondf$compounds
compoundorder<-cbind(compounds, c(1:length(compounds)))
data_replicates<-jsondf$replicates
data_outliers<-jsondf$outliers
col_num<-ncol(data_replicates)
row_num<-nrow(data_replicates)
rep_num<-dim(data_replicates)[3]
constraints<-jsondf$constrain
#Manualplate indicator, 0 by default, gets set to 1 if the data is a manual plate.
Manualplate<-0
#This function will test rearrange a manual plate's data to conform with the rest of the script.
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
#This blurb checks for whether or not a plate is a Manual plate and sets the information from before
#to correspond with that of a manual plate.
if(is.null(rep_num)==TRUE){
  data_replicates<-Manualdata(data_replicates)
  data_outliers<-Manualdata(data_outliers)
  rep_num<-dim(data_replicates)[3]
  col_num<-ncol(data_replicates)
  row_num<-nrow(data_replicates)
  Manualplate<-1
}
#Create empty array to check for outliers with. Array is either 3D or 2D, but will match data from JSON
data_replicates_checked<-array(NA, dim=dim(data_replicates))
#Check for outliers on the outlier grid. If it is an outlier, the data grid has the corresponding value
#set to NA
outliercheck<- function(replicates, outliers){
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
data_replicates_checked<-outliercheck(data_replicates, data_outliers)
#With outliers set as NA, bind the compound number to the first column in each plate and replace with numbered
#order 1-8. This differs for manual plates, so we check the indicator for a manual plate and decide
if(Manualplate == 0){
  data_replicates_cmpd<-abind(array(compoundorder[,2], replace(dim(data_replicates_checked),2,1)), 
  data_replicates_checked, along = 2)
}else if(Manualplate == 1){
  data_replicates_cmpd<-abind(array(rep(1:nrow(compoundorder), each = 3), dim=c(3,1,rep_num)),
  data_replicates_checked, along = 2)
  }
#Arrange the data into a single 2D array depending on 3D-ness of original array
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
#Convert the data columns into numeric and create average values, and place them into a table of the average values.
data_replicates[,2:ncol(data_replicates)]<-lapply(data_replicates[,2:ncol(data_replicates)], function(x) as.numeric(as.character(x)))
data_avgcounts<-data_replicates %>% group_by(V1)%>%summarise_all(mean, na.rm=TRUE)%>%
  arrange(match(V1,c(compoundorder[,2])), desc(V1))
#Prep assay information from JSON, most importantly the constant used in the equation. We term this constant
#as K-value. Someone with sufficient knowledge of the binding equation can tell you how and what it does. 
#Just know it as a constant that is used in the equation derived from radioactivity of the hot ligand.
curies<-as.numeric(jsondf$hot_count)/(2.22*10^12)
mmols<-curies/as.numeric(jsondf$hot_activity)
hotnM<-mmols*10^6/(as.numeric(jsondf$hot_volume)*10^-6)
Kd<-as.numeric(jsondf$dissociation_constant)
kvalue<-log10(1+(hotnM/Kd))
#Prep concentration/X-values from JSON. We keep the old X-values as we do need them all, but we drop
#any X-values that we can't use for regression. For instance, if a set of Y-values are all removed due to outliers,
#the corresponding X-value must be removed as well.
initialconcentrations<-as.data.frame(array(as.numeric(jsondf$concentrations), dim = dim(jsondf$concentrations)))
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
#Global Regression requires stacking Y values and X values into one long vector..
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
#Global regression will require indicator matrices to pair individual Ki's to 
#their respective data sets. However, for shared parameters this is not needed. BUT for shared parameters,
#any individual parameters must have their indicators removed so they don't conflict.
#We make a giant list where each vector is going to be a single indicator vector of 1's and 0's
#These vectors are similar to design matrices found in linear regression/ANOVA design matrices.
#1 indicates this data set value is for this parameter, while 0 indicates this value is not for this parameter.
Indicatorlist<-list()
for (i in 1:length(compounds)){
  Indicatorlist[[i]]<-rep(0, length(Yavg))
}
fullindicator<-rep(1, length(Yavg))
index<-0
for (i in 1:length(compounds)){
  indexprev<-index+1
  index<-index+(length(Yvals[[i]][2,])-1)
  Indicatorlist[[i]][indexprev:index]<-1
}
#Currently, the full indicator is only used for the only parameter that can be individualized, the "Bottom".
#Here we consider which Bottoms must be constrained, and create an individual parameter for those that are not
#constrained to be the same. Otherwise we use the Global Bottom.
constraintsBottoms<-as.list(compounds)
for (i in 1:length(constraints)){
  if (constraints[i] == 0){
    constraintsBottoms[[i]][["Bottom"]]<-noquote(paste0(quote(IndivBottom), i))
    constraintsBottoms[[i]][["Indicatorlist"]]<-noquote(paste0("Indicatorlist[[",i,"]]"))
    constraintsBottoms[[i]][["Indicator Number"]]<-i
    fullindicator<-fullindicator-Indicatorlist[[i]]
  } else if(constraints[i] == 1){
    constraintsBottoms[[i]][["Bottom"]]<-noquote(quote("GlobalBottom"))
    constraintsBottoms[[i]][["Indicatorlist"]]<-noquote(quote("fullindicator"))
    constraintsBottoms[[i]][["Indicator Number"]]<-0
  }
}
Bottomparams<-c()
for(i in 1:length(constraintsBottoms)){
  Bottomparams[i]<-constraintsBottoms[[i]][["Bottom"]]
}
Botindicatorsparams<-c()
for(i in 1:length(constraintsBottoms)){
  Botindicatorsparams[i]<-constraintsBottoms[[i]][["Indicatorlist"]]
}
Botindicatorsparams<-unique(Botindicatorsparams)
Botindicatorsparams<-sort(Botindicatorsparams)
Bottomparams<-unique(Bottomparams)
Bottomparams<-sort(Bottomparams)
#Now we construct the function. Basically we take all our parameters and put them into the equation 
#specified in the protocol book for Secondary Binding. We have to create this function programatically because
#the number of parameters differs across all data sets, so the function's body must be created according to how
#many parameters are available to be estimated.
BottomparamsxIndics<-c()
for(i in 1:length(Bottomparams)){
  BottomparamsxIndics[i]<-paste0("params$",Bottomparams[i],"*",Botindicatorsparams[i])
}
paramsxIndics<-c()
for (i in 1:length(compounds)){
  paramsxIndics[i]<-noquote(paste0("params$logKi",i,"*Indicatorlist[[",i,"]]"))
}
#These variables are expressions that are captured by the expr() function. !! is an indicator for unquoting.
s<-paste0("(",BottomparamsxIndics,collapse = "+",")")
t<-parse_expr((s))
r<-paste0(paramsxIndics,collapse = "+")
q<-parse_expr(r)
body<-expr((!!t)+((params$Top-(!!t))/(1+10^(xx-(!!q)-kvalue))))
#Generate Logki Parameter table. First set all ki's to -7. Then add top and bottom estimates for all tops and bottoms.
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
for (i in 1:length(compounds)){
  if (constraintsBottoms[[i]][["Indicator Number"]]!=0){
    paramslist[[paste0(quote(IndivBottom),i)]]<-bot_est
  }
  else if (constraintsBottoms[[i]][["Indicator Number"]] == 0){
    paramslist[["GlobalBottom"]]<-bot_est
  }
}

#Now we need to insert the correct number of variables into our formula. A huge pain. Earlier we created the 
#argument list, and body. new_function() is a function that creates a function from a given list of arguments
#and  a function body.
predbod<-body
predargs<-alist(params=1, xx=2)
g<-new_function(predargs, predbod)
#Create simulated values based on equation and starting estimates.
#simDnoisy <- g(paramslist, Xavg)
#The residual function used to gauage our residuals
residFun<- function(p, observed, xx) {
  observed - g(p,xx)
}
#The actual regression function. Literally one line. From the nls.lm package. It needs a list of parameters
#to be estimated, a residual function, a list of observations and corresponding X-values.
nls.out<- nls.lm(paramslist, fn = residFun, observed = Yavg, xx = Xavg)
#summary(nls.out)
#Compile results.
paramslisttotal<-list()
for(i in 1:length(compounds)){
  if (constraints[i] == 0){
    paramslisttotal[[i]]<-list(logKi=coef(nls.out)[[paste0("logKi",i)]],
                               Top=coef(nls.out)[["Top"]],Bottom=coef(nls.out)[[paste0("IndivBottom",i)]])
  } else if(constraints[i] == 1){
    paramslisttotal[[i]]<-list(logKi=coef(nls.out)[[paste0("logKi",i)]],
                               Top=coef(nls.out)[["Top"]],Bottom=coef(nls.out)[["GlobalBottom"]])
  }
}
#Grab predicted Y-values using the updated parameters of the equation.
getPredsingle<- function(logki, top, bottom, xx) {
  ((bottom) + ((top-bottom)/
  (1+(10^(xx-logki-kvalue)))))
}

output_assembly<-paramslisttotal
output<-list()
#Rename output_assembly using compound names. Now we can construct our JSON to be passed back.
names(output_assembly)<-c(compounds)
#Get predicted Y-values using initial X-values and new parameters. This is for the graph.
for(i in 1:length(compounds)){
output[[i]]<-getPredsingle(paramslisttotal[[i]][["logKi"]], paramslisttotal[[i]][["Top"]],
                           paramslisttotal[[i]][["Bottom"]], initialconcentrations[i,])
}
for(i in 1:length(compounds)){
  output_assembly[[i]][["Pred Y-Value"]]<-output[[i]]
}
#Name all of our output list elements
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
  finaloutput[[6]][[i]]<-as.vector(as.numeric(initialconcentrations[i,]))
  #finaloutput[[6]][[i]]<-as.vector(as.numeric(Xvals[[i]][,2:length(Xvals[[i]])]))
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

