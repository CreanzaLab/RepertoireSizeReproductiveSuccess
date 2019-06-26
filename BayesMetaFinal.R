########
#Coded by Cristina Robinson
#Last Modified 10-5-2018
#R Version 3.5.2
#JAGS v4
#metafor_v2.0-0    BEST_0.5.1   HDInterval_0.2.0
#phytools_0.6-44   dplyr_0.7.6  MCMCglmm_2.26
#maps_3.3.0        ape_5.1      coda_0.19-1
#rjags_4-6        xtable_1.8-3
########



#NOTE ON MULTREE!!!!
#this function will always say that the model failed to converge for our equations.  Do not trust it.
#Look at the effective sampling sizes.  The sample size for SE is 0, because it is fixed at 1, and this
#triggers mulTree to think that the model failed when everything is fine.  Trust the raw values!





##Set directory here
setwd("D:/Documents/R/2018-Meta-Analysis")
### Loads packages and my functions
source("BayesMetaFunctions.R")
dir <- getwd();Mainfolder <- "Output"
tag <- ""
#setwd(paste0(getwd(),"/Revision"))

######################
#get the consensus Tree and the multitree, fix names for analysis
#get the Field data, flip correlates on time so that large cor = success,
#move hatch to laying b/c only one hatch sample post-processing;
#variances are not appropriately sampled for N=1
set.seed(49)
MultiTree <- read.nexus("NewTree.nex")[sample(1000,100)]
TreeAinv <- GetRelatednessMatrix(file.path(dir, "NewTree.nex"))
RawData <- read.csv(file.path(dir, "AnalysisData.csv"), header = TRUE)
RawData[,"animal"] <- RawData$Species
flipsign <- which(RawData$MClass == "Time")
RawData$Raw.Cor[flipsign] <- RawData$Raw.Cor[flipsign]*-1
RawData$MType[which(RawData$MType == "Hatch")] <- "Laying"
RawData$MType <- droplevels(RawData$MType)

#set up the larger data sets
FullDataSet <- MakeSet(RawData, 38, "Full Dataset")
p.varfull <- var(FullDataSet$Zs)
OCDataSet <- MakeSet(RawData[which(is.na(RawData$O.C) == FALSE),], 38, "Song Stability Dataset")
p.var <- var(OCDataSet$Zs)
SylDataSet <- MakeSet(RawData[which(is.na(RawData$Syllable.Rep) == FALSE),], 38,"Syllable Repertoire Dataset")
p.var2 <- var(SylDataSet$Zs)


#make folder for plots
if(dir.exists(file.path(dir,"Output"))){
  if(file.exists(file.path(dir,"Output",tag))){
    setwd(file.path(dir, "Output",tag))
  }else{
    dir.create(file.path(dir, "Output",tag))
    setwd(file.path(dir, "Output",tag))
  } 
}else{
  dir.create(file.path(dir, "Output"))
  dir.create(file.path(dir, "Output",tag))
  setwd(file.path(dir, "Output",tag))
}



pdf("RepertoireSizeVersusZsasaskedforbyReviewer5or7dependingonhowyoucount.pdf")
#colors <- c('red', 'orange', 'yellow', 'green', 'blue', 'purple', 'black')
plot(log(SylDataSet$Syllable.Rep), SylDataSet$Zs,
     xlab="Log Average Species Syllable Repertoire Size",
     ylab="Fisher's Z", font.lab=2, pch=19,
     #col=colors[SylDataSet$MType],
     xaxt='n')
#legend('bottomright', levels(SylDataSet$MType), col=colors, pch=19)
abline(h=0)
axis(side=1, seq(5,1500, by=20), at=log(seq(5,1500, by=20)))
dev.off()



#FunnelPlot Series
pdf("Funnel Main.pdf");par(mfrow=c(2,2),mgp=c(1.75,.5,0), mar=c(4,4,.5,1))
FunnelPlots(FullDataSet, Groups = "Rep",1)
FunnelPlots(FullDataSet, Groups = "OC",2);dev.off()
#Test for Asymetry
sink("Async.txt")
print("Full");AsyncTests(FullDataSet$Zs, FullDataSet$Se)
print("Syl");AsyncTests(SylDataSet$Zs, SylDataSet$Se)
print("OC");AsyncTests(OCDataSet$Zs, OCDataSet$Se);sink(NULL)
#Distribution plot
pdf("Distribution of Syllable Repertoire Sizes.pdf");par(mfrow=c(1,1), mar=c(1,3,1,1), mgp=c(1.5,.5,0));DistributionPlotV2(FullDataSet);dev.off()






#Variance Testing, used the OC Data set, as this allows for testing all three Fixed effects models
#models run for 200K itt with burn in of 30K.  Variance estimates for prios was variances in the
#data/number of random effects tested (not incluing SE)



Threshs <- sort(unique(OCDataSet$Syllable.Rep))
Threshs <- Threshs[2:(length(Threshs)-1)]
FixedEffect <- c("Zs~1", "Zs~lnRep", "Zs~O.C-1", rep("Zs~SylRepBi-1", length(Threshs)))
Captions <- c(paste0(c("Population", "Song stability",
            rep("Repertoire size",length(Threshs))),
            " variance in the song stability dataset. ",
            "Labels the same as in Supplemental Table S1."))
Captions[4:length(Captions)] <- paste0(Captions[4:length(Captions)], "Threshold $>$=", Threshs, ".")
Threshs <- c(38,38, 38,Threshs)
OCDataSet["lnRep"] <- log(OCDataSet$Syllable.Rep)
for(i in seq_along(FixedEffect)){
  OCDataSet[,"SylRepBi"] <- cut(OCDataSet$Syllable.Rep, c(0,Threshs[i],2000), labels=c("Smaller", "Larger"), right=FALSE)

  set.seed(49)
  model1 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~Species+idh(Se):units", Vari = p.var)
  model2 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~MType+idh(Se):units", Vari = p.var)
  model3 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~Study+idh(Se):units", Vari = p.var)
  model4 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~Species+animal+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
  model5 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~MType+Study+idh(Se):units", Vari = p.var)
  model6 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~MType+Species+animal+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
  model7 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~Study+Species+animal+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
  model8 <- BayesMeta(Data=OCDataSet, Fix=eval(FixedEffect[i]), Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)

  VarianceTests <- list(model1, model2, model3,
                        model4, model5, model6,
                        model7, model8)
  QuickLatexTable(VarianceTable(VarianceTests, OCDataSet$Se, FALSE)[,2:4],
                  Captions[i], file.path(getwd(), "VarianceTables.txt"))
}
##Testing the metrics to see how they are different

set.seed(49)
model9 <- BayesMeta(Data=FullDataSet, Fix="Zs~MType-1", Ran="~Species+animal+Study+idh(Se):units", Vari = p.varfull, TreeA=TreeAinv)

QuickLatexTable(MetaTable(list(FullDataSet), list(model9), "MType", PSort=TRUE),
                 "Metrics model meta-analysis in the song stability dataset", file.path(getwd(), "MetricTable.txt"))









#Effects in the Full Population (Dataset is what changes)
set.seed(49)
model10.1 <- BayesMeta(Data=FullDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.varfull, TreeA=TreeAinv)
model10.2 <- BayesMeta(Data=SylDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
model10.3 <- BayesMeta(Data=OCDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)

makeChains(model10.1, Title="FullPop-FullData", Data=FullDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.varfull, TreeA=TreeAinv)
makeChains(model10.2, Title="FullPop-SylData", Data=SylDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
makeChains(model10.3, Title="FullPop-OCData", Data=OCDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)

ConvergePlot(model10.1, "FullPop-FullData")
ConvergePlot(model10.2, "FullPop-SylData")
ConvergePlot(model10.3, "FullPop-OCData")



QuickLatexTable(MetaTable(list(FullDataSet, SylDataSet, OCDataSet),
                          list(model10.1, model10.2, model10.3), "All", TRUE),
                 "Population meta-analysis", file.path(getwd(), "MainMetTable.txt"))



#with multiple different phylogenic trees
set.seed(49)
BayesMetaMultiTree(Data=FullDataSet, MultiTree=MultiTree, Name="Population", Fix="Zs~1",
                   Ran="~MType+Species+animal+Study+idh(Se):units",
                   Vari = p.varfull)



#Stability testing in stability dataset
set.seed(49)
model11 <- BayesMeta(Data=OCDataSet, Fix="Zs~O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)

makeChains(model11, Title="SongStab", Data=OCDataSet, Fix="Zs~O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)

ConvergePlot(model11, "SongStab")

QuickLatexTable(MetaTable(list(OCDataSet),
                          list(model11), "O.C", TRUE),
                "Song stability meta-analysis", file.path(getwd(), "MainMetTable.txt"))
QuickLatexTable(BESTPrinter(OCDataSet, "O.C"),
                "Song stability BEST Results in the song stability dataset", file.path(getwd(), "MainMetTable.txt"))
pdf("ForestOC.pdf", width=7, height=9.5);BayesForestPlot(OCDataSet, Group="OC", title="");dev.off()

set.seed(49)
Inc <- rep("Not Increasing", nrow(OCDataSet))
Inc[which(OCDataSet$Type %in% c("Increasing", "Increasing and also modifying", "Non-significant increases"))] <- "Increasing"
OCDataSet[,"Increasing"] <- Inc
model11.1 <- BayesMeta(Data=OCDataSet, Fix="Zs~Increasing-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
#ConvergePlot(model11.1, "Increasing1")
QuickLatexTable(MetaTable(list(OCDataSet), list(model11.1), "Increasing", TRUE),
                "Song stability increasing less stringent meta-analysis", file.path(getwd(), "NewContResults.txt"))
QuickLatexTable(BESTPrinter(OCDataSet, "Increasing"),
                "Song stability BEST Results in the song stability dataset when only species that are known to increase their repertoires are considered plastic.", file.path(getwd(), "NewContResults.txt"))

Inc <- rep("Not Increasing", nrow(OCDataSet))
Inc[which(OCDataSet$Type %in% c("Increasing", "Increasing and also modifying"))] <- "Increasing"
OCDataSet[,"Increasing"] <- Inc
model11.2 <- BayesMeta(Data=OCDataSet, Fix="Zs~Increasing-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
#ConvergePlot(model11.2, "Increasing2")

QuickLatexTable(MetaTable(list(OCDataSet), list(model11.2), "Increasing", TRUE),
                "Song stability increasing stringent meta-analysis", file.path(getwd(), "NewContResults.txt"))
QuickLatexTable(BESTPrinter(OCDataSet, "Increasing"),
                "Song stability BEST Results in the song stability dataset when only species that are known to increase their repertoires are considered plastic.", file.path(getwd(), "NewContResults.txt"))

#with multiple different phylogenic trees
set.seed(49)
BayesMetaMultiTree(Data=OCDataSet, MultiTree=MultiTree, Name="SongStab", Fix="Zs~O.C-1",
          Ran="~MType+Species+animal+Study+idh(Se):units",
          Vari = p.var)





#All Thresholds in the syllable repertoire datset, Seed set internally
ThresholdTest <- JackThresh(SylDataSet, Vari=p.var2, Tree=TreeAinv, Type="Thresh")
Thresholds <- sort(unique(SylDataSet$Syllable.Rep))[2:(length(unique(SylDataSet$Syllable.Rep))-1)]
ThreshDataSets <- list()
for(i in seq_along(ThresholdTest)){
  #RHat done in JackThresh
  ConvergePlot(ThresholdTest[[i]], paste("SylRep", Thresholds[i]))
  SylDataSet[,"SylRepBi"] <- cut(SylDataSet$Syllable.Rep, c(0,Thresholds[i],2000), labels=c("Smaller", "Larger"), right=FALSE)
  ThreshDataSets[[i]] <- SylDataSet
}
QuickLatexTable(MetaTable(ThreshDataSets, ThresholdTest, "SylRepBi", Thresher=Thresholds),
                "Syllable repertoire model meta-analysis in the syllable repertoire dataset", file.path(getwd(), "MainMetTable.txt"))
QuickLatexTable(BESTPrinter(SylDataSet, "SylRepBi", Thresholds[2:(length(Thresholds))]),
            "Syllable repertoire BEST analysis in the syllable repertoire dataset", file.path(getwd(), "MainMetTable.txt"))
pdf("ForestREP.pdf", width=7, height=11);BayesForestPlot(SylDataSet, Group="Rep", title="");dev.off()





#Subset of syllable thresholds in OC dataset
Picks <- c(18.5, 38, 216)
set.seed(49)
DataSets <- list();ModelList <- list()
for(i in seq_along(Picks)){
  OCDataSet$SylRepBi <-  cut(OCDataSet$Syllable.Rep, c(0,Picks[i],2000), labels=c("Smaller", "Larger"), right=FALSE) 
  DataSets[[i]] <- OCDataSet
  ModelList[[i]] <- BayesMeta(Data=DataSets[[i]], Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
  
  makeChains(ModelList[[i]], Title=paste("Syl in song", Picks[i]), Data=DataSets[[i]], Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
             
  ConvergePlot(ModelList[[i]], paste("SylRepOCDataset", Picks[i], sep="-"))
}
QuickLatexTable(MetaTable(DataSets,ModelList, "SylRepBi", Thresher = Picks),
                "Syllable repertoire model meta-analysis in the Song Stability dataset", file.path(getwd(), "SupMetaTables.txt"))


for(i in seq_along(Picks)) {
  #with multiple different phylogenic trees
  set.seed(49)
  SylDataSet$SylRepBi <-  cut(SylDataSet$Syllable.Rep, c(0, Picks[i], 2000), labels=c("Smaller", "Larger"), right=FALSE)
  BayesMetaMultiTree(Data=SylDataSet, MultiTree=MultiTree, Name=paste("SylDis", Picks[i]), Fix="Zs~SylRepBi-1",
                     Ran="~MType+Species+animal+Study+idh(Se):units",
                     Vari=p.var2)
}




#Everything with RepSize Continuous
set.seed(49)
SylDataSet["lnRep"] <- log(SylDataSet$Syllable.Rep)
model13 <- BayesMeta(Data=SylDataSet, Fix="Zs~lnRep", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv, PR=TRUE)

makeChains(model13, Title="continuous", Data=SylDataSet, Fix="Zs~lnRep", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv, PR=TRUE)

ConvergePlot(model13, "Continuous model")
QuickLatexTable(MetaTable(list(SylDataSet), list(model13), "All", Groups=c("Intercept", "Larger")),
                "Continuous syllable repertoire model meta-analysis in the repertoire dataset", file.path(getwd(), "MainMetTable.txt"))

#continuous species jackknife
set.seed(49)
Spec <- unique(SylDataSet$Species)
DataSets <- vector(length = length(Spec), mode = "list")
Models <- vector(length = length(Spec), mode = "list")
for(i in seq_along(Spec)){
  DataSets[[i]] <- SylDataSet[-which(SylDataSet$Species == Spec[i]),]
  tempvar <- var(DataSets[[i]]$Zs)
  Models[[i]] <- BayesMeta(Data=DataSets[[i]], Fix="Zs~lnRep", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = tempvar, TreeA=TreeAinv)
}
ContSpecJack <- MetaTable(DataSets, Models, "All", Groups=c("Intercept", "Larger"))
ContSpecJack["Removed"] <- rep(Spec, each=2,1)
ContSpecJack <- ContSpecJack[,c(7, 1:6)]

QuickLatexTable(ContSpecJack,
                paste0("Continuous syllable repertoire model meta-analysis with each species removed in turn."),
                file.path(getwd(), "NewContResults.txt"))
#continuous, remove top and bottoms
set.seed(49)
RemoveExtremes <- SylDataSet[order(SylDataSet$Syllable.Rep),]
Spec2 <- unique(RemoveExtremes$Species)
BottomData <- vector(length=4, mode="list")
TopData <- vector(length=4, mode="list")
n <- 3:9
for(i in n){
  Tops <- which(RemoveExtremes$Species %in% Spec2[(length(Spec2)-i):length(Spec2)])
  Bottom <- which(RemoveExtremes$Species %in% Spec2[1:i])
  VarTop <- var(RemoveExtremes$Zs[-Tops])
  VarBottom <- var(RemoveExtremes$Zs[-Bottom])
  model13Top <- BayesMeta(Data=RemoveExtremes[-Tops,], Fix="Zs~lnRep", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = VarTop, TreeA=TreeAinv)
  model13Bottom <- BayesMeta(Data=RemoveExtremes[-Bottom,], Fix="Zs~lnRep", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = VarBottom, TreeA=TreeAinv)
  TopData[[i-2]] <- MetaTable(list(RemoveExtremes[-Tops,]),
                                 list(model13Top), "All", Groups=c("Intercept", "Larger"))
  BottomData[[i-2]] <- MetaTable(list(RemoveExtremes[-Bottom,]),
                                list(model13Bottom), "All", Groups=c("Intercept", "Larger"))
}
TopData <- do.call(rbind.data.frame, TopData)
BottomData <- do.call(rbind.data.frame, BottomData)
TopData["# Species Removed"] <- rep(n, each=2,1)
BottomData["# Species Removed"] <- rep(n, each=2,1)
TopData <- TopData[,c(7, 1:6)]
BottomData <- BottomData[,c(7, 1:6)]

QuickLatexTable(TopData,
                paste0("Continuous syllable repertoire model meta-analysis with # species with the largest syllable repertoires removed."),
                file.path(getwd(), "NewContResults.txt"))
QuickLatexTable(BottomData,
                paste0("Continuous syllable repertoire model meta-analysis with # species with the smallest syllable repertoires removed."),
                file.path(getwd(), "NewContResults.txt"))


PredictionData <- predict(model13,marginal=NULL, posterior="all",
                          type="response", interval = "confidence")
pdf("Prediction Data.pdf")
plot(x=SylDataSet$Zs, y=PredictionData[,"fit"], pch=19,
     xlab="Real Fisher's Z", ylab="Predicted Fisher's Z",
     xlim=c(-1,1.5), ylim=c(-1,1.5),
     panel.first={
       abline(a=0, b=1, lwd=2, lty=2, col="grey40")
       segments(x0=SylDataSet$Zs, x1=SylDataSet$Zs,
                y0=PredictionData[,"lwr"], y1=PredictionData[,"upr"], col="grey80")
       })
dev.off()



#with multiple different phylogenic trees
set.seed(49)
BayesMetaMultiTree(Data=SylDataSet, MultiTree=MultiTree, Name="SylCont", Fix="Zs~lnRep",
                   Ran="~MType+Species+animal+Study+idh(Se):units",
                   Vari = p.var2)





#simulate
set.seed(NULL)
for (i in seq_along(1:16)){
  pred_means <- vector()
  predictions <- simulate.MCMCglmm(glmm[[i]], nsim = 5000)
  #print(predictions)
  pred_means <- colMeans(predictions)
  #print(pred_means)
  #plot histogram of mean of predictions and line of mean of actual data
  print(hist(pred_means, main=glmm[[i]]$Fixed[[1]][[2]]))
  print(abline(v=mean(yearTimeRec_na_df[[i+10]])))
}






#SOng Stability and Repertoire size together
DataSets <- list();Models <- list();Models2 <- list()

for(i in seq_along(Picks)){
  OCDataSet[,"SylRepBi"] <- cut(OCDataSet$Syllable.Rep, c(0,Picks[i],2000), labels=c("Smaller", "Larger"), right=FALSE)
  OCDataSet[,"SylOc"] <- SetGroups(OCDataSet)
  DataSets[[i]] <- OCDataSet
  #set.seed(49)
  #Models[[i]] <- BayesMeta(Data=OCDataSet, Fix="Zs~SylRepBi:O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
  #makeChains(Models[[i]], Title=paste("Interaction",Picks[i]),Data=OCDataSet, Fix="Zs~SylRepBi:O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
  
  set.seed(49)
  Models2[[i]] <- summary(BayesMeta(Data=OCDataSet, Fix="Zs~SylRepBi+O.C", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv))
  makeChains(Models2[[i]], Title=paste("No Interaction",Picks[i]),Data=OCDataSet, Fix="Zs~SylRepBi:O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
  #ConvergePlot(Models[[i]], paste("RepandSongStab",Picks[i],sep="-"))
}
Groupss <- c("SmallerStable", "LargerStable", "SmallerPlastic", "LargerPlastic")
QuickLatexTable(MetaTable(DataSets, Models, "SylOc", Groups=Groupss, Thresher = rep(Picks, each=2)),
                "Combined syllable repertoire and song stability model meta-analysis in the song stability dataset",
                file.path(getwd(), "MainMetTable.txt"))







###new datasets Without Offspring and minus offspring and EPP
NoOffDataSet <- SylDataSet[which(SylDataSet$MClass != "Offspring"),]
NoOfforEPPDataSet <- NoOffDataSet[which(NoOffDataSet$MType != "EPP"),]
p.var3 <- var(NoOffDataSet$Zs)
p.var4 <- var(NoOfforEPPDataSet$Zs)


DataSetNoOff <- list();ModelsNoOff <- list()
DataSetNoOffEPP <- list();ModelsNoOffEPP <- list()
for(i in seq_along(Picks)){
  NoOffDataSet[,"SylRepBi"] <- cut(NoOffDataSet$Syllable.Rep, c(0,Picks[i],2000), labels=c("Smaller", "Larger"), right=FALSE)
  NoOfforEPPDataSet[,"SylRepBi"] <- cut(NoOfforEPPDataSet$Syllable.Rep, c(0,Picks[i],2000), labels=c("Smaller", "Larger"), right=FALSE)
  DataSetNoOff[[i]] <- NoOffDataSet 
  DataSetNoOffEPP[[i]] <- NoOfforEPPDataSet
  
  set.seed(49)
  ModelsNoOff[[i]] <- BayesMeta(Data=NoOffDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var3, TreeA=TreeAinv)
  ModelsNoOffEPP[[i]] <- BayesMeta(Data=NoOfforEPPDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var4, TreeA=TreeAinv)
  
  #ConvergePlot(ModelsNoOff[[i]], paste("NoOff", Picks[i], sep="-"))
  #ConvergePlot(ModelsNoOffEPP[[i]], paste("NoOfforEPP", Picks[i], sep="-"))
}
QuickLatexTable(MetaTable(DataSetNoOff, ModelsNoOff, "SylRepBi", Thresher = Picks),
                "Repertoire model meta-analysis in the dataset without Offspring Measurements",
                file.path(getwd(), "NoOffMetaTables.txt"))
QuickLatexTable(MetaTable(DataSetNoOffEPP, ModelsNoOffEPP, "SylRepBi", Thresher = Picks),
                "Repertoire model meta-analysis in the dataset without offspring or EPP Measurements",
                file.path(getwd(), "NoOffMetaTables.txt"))
QuickLatexTable(BESTPrinter(NoOffDataSet, "SylRepBi", Picks),
                "Syllable repertoire BEST analysis in the no offspring dataset", file.path(getwd(), "NoOffMetaTables.txt"))
QuickLatexTable(BESTPrinter(NoOfforEPPDataSet, "SylRepBi", Picks),
                "Repertoire model meta-analysis in the dataset without offspring or EPP Measurements", file.path(getwd(), "NoOffMetaTables.txt"))


SpeciesDatasets <- data.frame(Full=as.character(unique(FullDataSet$Species)),
                              Repertoire=rep("",27),
                              Stability=rep("",27),
                              NoOffspring=rep("",27),
                              NoOffspringorEPP=rep("",27),
                              stringsAsFactors=FALSE)
Syl <- as.character(unique(SylDataSet$Species))
SpeciesDatasets[which(SpeciesDatasets[,1] %in% Syl),2] <- Syl
OC <- as.character(unique(OCDataSet$Species))
SpeciesDatasets[which(SpeciesDatasets[,1] %in% OC),3] <- OC
NoOFF <-as.character(unique(NoOffDataSet$Species))
SpeciesDatasets[which(SpeciesDatasets[,1] %in% NoOFF),4] <- NoOFF
NoOffEPP <-as.character(unique(NoOfforEPPDataSet$Species))
SpeciesDatasets[which(SpeciesDatasets[,1] %in% NoOffEPP),5] <- NoOffEPP


print(xtable(SpeciesDatasets),include.rownames=FALSE)





###############
#Backupchecks

#Each of these has to run many meta-analysis and therefore take a long time!!!!

#Species jackknife, seed set internally
#Seed set internally
Species <- as.character(unique(SylDataSet$Species))
for(i in seq_along(Picks)){
  SylDataSet[,"SylRepBi"] <- cut(SylDataSet$Syllable.Rep, c(0,Picks[i],2000), labels=c("Smaller", "Larger"), right=FALSE)
  SPJack <- JackThresh(SylDataSet, Vari=p.var2, Tree=TreeAinv, Type="SpJack",Group = "Rep")
  SPDatasets <- list()
  for(j in seq_along(Species)){
    SPDatasets[[j]] <- SylDataSet[-which(SylDataSet$Species==Species[j]),]  
  }
  SPJackTable <- MetaTable(SPDatasets,SPJack, "SylRepBi")
  SPJackTable["Removed"] <- rep(Species, each=2,1)
  SPJackTable <- SPJackTable[,c(7, 1:6)]
  QuickLatexTable(SPJackTable,
                  paste0("Syllable repertoire model meta-analysis with each species removed in turn, threshold $>$=",Picks[i]),
                  file.path(getwd(), "JackKnifeTables.txt"))
}


#switching the plasticity categorization, seed set internally
OCJack <- JackThresh(OCDataSet, Vari=p.var, Tree=TreeAinv, Type="OCJack",Group = "OC")
Species <- as.character(unique(OCDataSet$Species))
OCDatasets <- list()
for(j in seq_along(Species)){
  OCCopy <- OCDataSet
  ind <- which(OCDataSet$Species==Species[j])
  if(OCDataSet$O.C[ind[1]] == "Stable"){
    OCCopy$O.C[ind] <- "Plastic"
  }else{
    OCCopy$O.C[ind] <- "Stable"
  }
  OCDatasets[[j]] <- OCCopy   
}
OCJackTable <- MetaTable(OCDatasets, OCJack, "O.C")
OCJackTable["Switched"] <- rep(Species, each=2,1)
OCJackTable <- OCJackTable[,c(7, 1:6)]
QuickLatexTable(OCJackTable,
                "Song stability model meta-analysis in song stability dataset with song stability categories of individual species switched",
                file.path(getwd(), "JackKnifeTables.txt"))





#not totally sure shat reviewer 3 wants
Smaller <- OCDataSet[which(OCDataSet$Syllable.Rep < 38),]
p.varsmall <- var(Smaller$Zs)
Larger <- OCDataSet[which(OCDataSet$Syllable.Rep >= 38),]
p.varlarge <- var(Larger$Zs)
set.seed(49)
Smalltest <- BayesMeta(Smaller, Fix="Zs~O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.varsmall, TreeA=TreeAinv)
Largetest <- BayesMeta(Larger, Fix="Zs~O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.varlarge, TreeA=TreeAinv)
summary(Smalltest)
summary(Largetest)

#Territory Control
#make folder for plots
if(file.exists(paste0(dir,"/Output/",tag, "/Terr"))){
  setwd(file.path(dir, "Output",tag, "Terr"))
} else {
  dir.create(file.path(dir, "Output",tag, "Terr"))
  setwd(file.path(dir, "Output",tag, "Terr"))
  
}

#Chnage out the dataset
Terr <- which(is.na(RawData$Territory.Control) == FALSE)
RawTerr <- RawData
RawTerr$Raw.Cor[Terr] <- RawTerr$Territory.Control[Terr]
Flip <- which(RawTerr$MClass[Terr] == "Time")
RawTerr$Raw.Cor[Terr][Flip] <- RawTerr$Raw.Cor[Terr][Flip]*-1

#set up the larger data set
TFullDataSet <- MakeSet(RawTerr, 38, "Full Dataset")
p.varfullT <- var(TFullDataSet$Zs)
TOCDataSet <- MakeSet(RawTerr[which(is.na(RawData$O.C) == FALSE),], CUT=38, "Song Stability Dataset")
p.var5 <- var(TOCDataSet$Zs)
TSylDataSet <- MakeSet(RawTerr[which(is.na(RawData$Syllable.Rep) == FALSE),], 38, "Syllable Repertoire Dataset")
p.var6 <- var(TSylDataSet$Zs)



#FunnelPlot Series
pdf("Funnel Territory.pdf");par(mfrow=c(2,2),mgp=c(1.75,.5,0), mar=c(4,4,.5,1))
FunnelPlots(TFullDataSet, Groups = "Rep",1)
FunnelPlots(TFullDataSet, Groups = "OC",2);dev.off()

#Test for Asymetry
sink("TAsync.txt")
print("Full");AsyncTests(TFullDataSet$Zs, TFullDataSet$Se)
print("Syl");AsyncTests(TSylDataSet$Zs, TSylDataSet$Se)
print("OC");AsyncTests(TOCDataSet$Zs, TOCDataSet$Se);sink(NULL)








#POPULATION
set.seed(49)
model1.1T <- BayesMeta(Data=TFullDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.varfullT, TreeA=TreeAinv)
model1.2T <- BayesMeta(Data=TSylDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var6, TreeA=TreeAinv)
model1.3T <- BayesMeta(Data=TOCDataSet, Fix="Zs~1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var5, TreeA=TreeAinv)
ConvergePlot(model1.1T, "PopTFull")
ConvergePlot(model1.2T, "PopTSyl")
ConvergePlot(model1.3T, "PopTOC")
QuickLatexTable(MetaTable(list(TFullDataSet, TSylDataSet, TOCDataSet),
                          list(model1.1T, model1.2T, model1.3T), "All"),
                "Population Meta-Analysis with Territory-controlled Measurements",
                file.path(getwd(),"TablesT.txt"))








#RepSize Dataset
set.seed(49)
TDataSets <- list(); TModels <- list()
for(i in seq_along(Picks)){
  TSylDataSet[,"SylRepBi"] <- cut(TSylDataSet$Syllable.Rep, c(0,Picks[i],2000), labels=c("Smaller", "Larger"), right=FALSE)
  TDataSets[[i]] <- TSylDataSet
  TModels[[i]] <- BayesMeta(Data=TSylDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var6, TreeA=TreeAinv)
  ConvergePlot(TModels[[i]],paste("RepSizeT", Picks[i],sep="-"))
}
QuickLatexTable(MetaTable(TDataSets, TModels, "SylRepBi", Thresher = Picks),
                "Syllable repertoire model in the territory-controlled repertoire dataset meta-analysis",
                file.path(getwd(),"TablesT.txt"))
QuickLatexTable(BESTPrinter(TSylDataSet, "SylRepBi", Picks),
                "Syllable repertoire BEST results with territory-controlled measurements", file.path(getwd(),"TablesT.txt"))
pdf("ForestREPT.pdf", width=7, height=11);BayesForestPlot(TSylDataSet, Group="Rep", title="");dev.off()


set.seed(49)
TSylDataSet["lnRep"] <- log(TSylDataSet$Syllable.Rep)
model13T <- BayesMeta(Data=TSylDataSet, Fix="Zs~lnRep", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var6, TreeA=TreeAinv)
ConvergePlot(model13T, "Continuous model")
QuickLatexTable(MetaTable(list(SylDataSet), list(model13), "All", Groups=c("Intercept", "Larger")),
                "Continuous syllable repertoire model meta-analysis in the repertoire dataset with territory-controlled measurements.",
                file.path(getwd(), "NewContResults.txt"))




#Stability testing
set.seed(49)
model3T <- BayesMeta(Data=TOCDataSet, Fix="Zs~O.C-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var5, TreeA=TreeAinv)

ConvergePlot(model3T, "SongStabT")
QuickLatexTable(MetaTable(list(TOCDataSet),list(model3T),"O.C"),
                "Song stability model in the territory-controlled song stability dataset meta-analysis", file.path(getwd(),"TablesT.txt"))
QuickLatexTable(BESTPrinter(TOCDataSet, "O.C", Groups=c("Stable", "Plastic")),
                            "Song stability BEST results with territory-controlled measurements", file.path(getwd(),"TablesT.txt"))
pdf("ForestOCT.pdf", width=7, height=9.5);BayesForestPlot(TOCDataSet, Group="OC", title="");dev.off()








#make folder for plots
if(file.exists(paste0(dir,"/Output/",tag, "/Misc"))){
  setwd(file.path(dir, "Output",tag, "Misc"))
} else {
  dir.create(file.path(dir, "Output",tag, "Misc"))
  setwd(file.path(dir, "Output",tag, "Misc"))
  
}

###Plots and analysis for MinMax
pdf("Distribution of MinMax Syllable Repertoire Sizes.pdf");par(mfrow=c(1,2))
MinMaxDisPlot(SylDataSet);dev.off()


SylDataSetMAX <- SylDataSet
SylDataSetMAX$Syllable.Rep <- SylDataSet$Syl.Max
ThresholdTestMAX <- JackThresh(SylDataSetMAX, Vari=p.var2, Tree=TreeAinv, Type="Thresh")
ThresholdsMAX <- sort(unique(SylDataSetMAX$Syllable.Rep))[2:(length(unique(SylDataSetMAX$Syllable.Rep))-1)]
MaxDataSets <- list()
for(i in seq_along(ThresholdsMAX)){
  SylDataSetMAX[,"SylRepBi"] <- cut(SylDataSetMAX$Syllable.Rep, c(0,ThresholdsMAX[i],2000), labels=c("Smaller", "Larger"), right=FALSE)
  MaxDataSets[[i]] <- SylDataSetMAX
}

SylDataSetMIN <- SylDataSet
SylDataSetMIN$Syllable.Rep <- SylDataSet$Syl.Min
ThresholdTestMIN <- JackThresh(SylDataSetMIN, Vari=p.var2, Tree=TreeAinv, Type="Thresh")
ThresholdsMIN <- sort(unique(SylDataSetMIN$Syllable.Rep))[2:(length(unique(SylDataSetMIN$Syllable.Rep))-1)]
MinDataSets <- list()
for(i in seq_along(ThresholdsMIN)){
  SylDataSetMIN[,"SylRepBi"] <- cut(SylDataSetMIN$Syllable.Rep, c(0,ThresholdsMIN[i],2000), labels=c("Smaller", "Larger"), right=FALSE)
  MinDataSets[[i]] <- SylDataSetMIN
}

QuickLatexTable(MetaTable(MaxDataSets,ThresholdTestMAX, "SylRepBi", Thresher = ThresholdsMAX),
                "Syllable repertoire model meta-analysis in the max syllable repertoire dataset",
                file.path(getwd(), "MinMax.txt"))
QuickLatexTable(MetaTable(MinDataSets,ThresholdTestMIN, "SylRepBi", Thresher = ThresholdsMIN),
                "Syllable repertoire model meta-analysis in the min syllable repertoire dataset",
                file.path(getwd(), "MinMax.txt"))
QuickLatexTable(BESTPrinter(SylDataSetMAX, "SylRepBi", ThresholdsMAX[2:length(ThresholdsMAX)]),
                "Syllable repertoire BEST analysis in the max syllable repertoire dataset",
                file.path(getwd(), "MinMax.txt"))
QuickLatexTable(BESTPrinter(SylDataSetMIN, "SylRepBi", ThresholdsMIN[2:length(ThresholdsMIN)]),
                "Syllable repertoire BEST analysis in the min syllable repertoire dataset",
                file.path(getwd(), "MinMax.txt"))


set.seed(49)
SylDataSetMAX$lnRep <- log(SylDataSetMAX$Syllable.Rep)
SylDataSetMIN$lnRep <- log(SylDataSetMIN$Syllable.Rep)
model13MAX <- BayesMeta(Data=SylDataSetMAX, Fix="Zs~lnRep", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
model13MIN <- BayesMeta(Data=SylDataSetMIN, Fix="Zs~lnRep", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)

MinMaxCont <- MetaTable(list(SylDataSetMAX, SylDataSetMIN), list(model13MAX,model13MIN), "All", Groups=c("Intercept", "Larger"))
MinMaxCont["Min or Max"] <- rep(c("Max", "Min"), each=2)
MinMaxCont <- MinMaxCont[,c(7,1:6)]
QuickLatexTable(MinMaxCont,
                "Continuous syllable repertoire model meta-analysis in the max and min syllable repertoire datasets.",
                file.path(getwd(), "NewContResults.txt"))

#Just the lowest versus the highest
EndsDataSet <- SylDataSet
Remove <- which((EndsDataSet$Syllable.Rep > 20 & EndsDataSet$Syllable.Rep < 100))
EndsDataSet <- EndsDataSet[-Remove,]
set.seed(49)
EndsModel <- BayesMeta(Data=EndsDataSet, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)

QuickLatexTable(MetaTable(list(EndsDataSet), list(EndsModel), "SylRepBi"),
                "Syllable repertoire BEST analysis in the ends only dataset",
                file.path(getwd(),"JustforOriginalRev2.txt"))
QuickLatexTable(BESTPrinter(EndsDataSet, "SylRepBi",100),
                "Syllable repertoire BEST analysis in the min syllable repertoire dataset",
                  file.path(getwd(),"JustforOriginalRev2.txt"))




#Graph of song and SYl rep correlation
Table <- RawData[which(!duplicated(RawData$Species)),]
remove <- which(is.na(Table$Song.Rep))
Data <- data.frame(x=log(Table$Syllable.Rep[-remove]), y=log(Table$Song.Rep[-remove]))
pdf("SylSongRep.pdf")
par(mgp=c(3,.5,0), mar=c(4,4,1,1), mfrow=c(2,2))
plot(Data, xaxt="n", yaxt="n",
     xlab="Syllable Repertoire Size", ylab="Song Repertoire Size", pch=20,
     font.lab=2)
axis(1, 0:7, round(exp(0:7), digits=1), las=2)
axis(2, 0:7, round(exp(0:7), digits=1), las=1)
Results <- cor.test(Data$x, Data$y, alternative ="greater", method="spearman")
text(x=4,y=.5, adj=0,
     paste0("Rho = ", round(Results$estimate, digits=3),
            "\np-value = ", round(Results$p.value, digits=3)), cex=.8)
dev.off()



#ARBITRARY DATA
#Generating this data takes ~7.5 hours.
######
#Generate the data
Nitt <- 500
ArbitraryDataSets <- GenerateArbitraryDataSets(SylDataSet,Nitt)
ArbitraryModels <- list()
for(i in 1:Nitt){
  set.seed(49)
  ArbitraryModels[[i]] <- BayesMeta(Data=ArbitraryDataSets[[i]], Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var2, TreeA=TreeAinv)
}
ArbitraryBESTResults <- data.frame(character(Nitt),character(Nitt),character(Nitt), stringsAsFactors = FALSE)
for(i in 1:Nitt){
  ArbitraryBESTResults[i,] <- as.vector(BESTPrinter(ArbitraryDataSets[[i]], "SylRepBi"))
}
colnames(ArbitraryBESTResults) <- c("BEST Mean", "95% CredInt", "%<0")

######
#Data reoganization; go here post loading rDATA
pdf("Arbitrary Runs.pdf");ArbitraryPlots(ArbitraryDataSets, ArbitraryModels, ArbitraryBESTResults, ThresholdTest, Thresholds);dev.off()


#ArbTable[,"#Species"] <- c(rbind(SmallSpec, length(Spec)-SmallSpec))
#ArbTable[,"#Measure"] <- c(rbind(SmallMeasure, nrow(Arbitrary)-SmallMeasure))
#ArbTable <- ForestTable(ArbModel, rep(list(Arbitrary), Nu), Picks=20)[5:7]
#SmallMeasure <- colSums(ArbClasses)
#SmallSpec <- sapply(1:Nu, function(x) length(unique(Arbitrary$Species[as.logical(ArbClasses[,x])])))
#ArbTable["#Species"] <- c(rbind(SmallSpec, length(Spec)-SmallSpec))
#ArbTable["#Measure"] <- c(rbind(SmallMeasure, nrow(Arbitrary)-SmallMeasure))
#ArbTable["Threshold"] <- NULL
#percentsig <- which(percent < 2.5 | percent > 97.5)
#sigmodels <- sig[which(sig %in% percentsig)]
#TabledBestResults <- matrix(unlist(BestResults), nrow=Nu, byrow = TRUE)[,2:4]





#which species in which group
#Models <- lapply(sigmodels, function(x) ArbModel[[x]])
#ArbTableSig <- ForestTable(Models, rep(list(Arbitrary), length(Models)), Picks=20)[,c(2,5:7)]
#Sigs <- which(ArbTableSig$pMCMC < .025)
#LSig <- which(Sigs%%2==0)

par(mar=c(11,3.5,0,1),mgp=c(1.5,.5,0))
hist(rep(1:length(Spec),SplitLarge), breaks=0:25, col=rgb(.5,.5,1,.8),
     ylim=c(0, 50), xaxt="n", xlab="", main="", yaxt="n",
     ylab="Percentage of Significant Groups With a Species        ",
     font.lab=2, cex.lab=.9)
rect(0,16.5,length(Spec),27.5, col=rgb(.4,.4,.4,.4), border=NA)
axis(1, seq(.5,length(Spec),by=1), paste0(SylOrder, "(",
                                          sort(Arbitrary[which(duplicated(as.character(Arbitrary$Species))==FALSE),]$Syllable.Rep),
                                          ")"), las=2, pos=0, cex.axis=.6, font.axis=3)
axis(2, seq(0,44,by=44/10), seq(0,100,by=10), las=2, pos=0)
mtext("B", side=2, las=2, at=44.5, line=1.5, font=2)
segments(0,44,25,44, lwd=1, lty=1)
segments(25,0,25,44, lwd=1, lty=1)
segments(0,22,25,22, lwd=2, lty=2)

