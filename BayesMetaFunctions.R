####
#Coded by Cristina Robinson
#Last Modified 10-5-2018
#R Version 3.5.2
#JAGS v4
#metafor_v2.0-0    BEST_0.5.1   HDInterval_0.2.0
#phytools_0.6-44   dplyr_0.7.6  MCMCglmm_2.26
#maps_3.3.0        ape_5.1      coda_0.19-1
#rjags_4-6        xtable_1.8-3
####


#NOTE ON MULTREE!!!!
#this function will always say that the model failed to converge for our equations.  Do not trust it.
#Look at the effective sampling sizes.  The sample size for SE is 0, because it is fixed at 1, and this
#triggers multree to think that the model failed when everything is fine.  Trust the raw values!

#The multree print outs are also slightly broken, because, again, our effective sample size for SE is 0.
#The info for the fixed effects and phylogenetic variance are correct, but the infor for residual varaance is wrong.

#load packages
Packs <- c("Matrix", "coda", "ape", "maps","MCMCglmm",
           "dplyr", "metafor", "xtable", "phytools",
           "rjags","BEST", "devtools", "caper")

for(i in Packs){
  if(!eval(parse(text=paste0("require(",i,")")))){
    install.packages(i)
  }
  eval(parse(text=paste0("library(",i,")")))
}

if(!require("multree")){
  install_github("TGuillerme/mulTree", ref = "release")
}

library("mulTree")

rm(list=objects())


#convert to Zs, set  datastructure and tree
GetRelatednessMatrix <- function(path){
  Tree <- read.nexus(path)
  Tree <- consensus.edges(Tree)
  #this is for finicky jerkface trees
  Tree <- multi2di(Tree, .0001)
  Tree <- force.ultrametric(Tree)
  Tree$edge.length[which(Tree$edge.length == 0)] <- .0001
  Tree <- force.ultrametric(Tree)
  Tree$tip.label <- gsub("_", " ", Tree$tip.label)
  return(inverseA(Tree, nodes="TIPS")$Ainv)
}
fischerTransform <- function(Rs, Ns){
  #correction for rtoz positive bias
  Rs <- Rs-((Rs*(1-Rs^2))/(2*(Ns-3)))
  Zs <- transf.rtoz(Rs)
  return(Zs)
}

MakeSet <- function(Data, CUT=38, Name=""){
  Data[,"SylRepBi"] <- cut(Data$Syllable.Rep, c(0,CUT,2000), labels=c("Smaller", "Larger"), right=FALSE)
  Data[,"Zs"] <- fischerTransform(Data$Raw.Cor, Data$nBirds)
  Data[,"Se"] <- 1/sqrt(Data$nBirds-3)
  Data[,"All"] <- rep(Name, nrow(Data))
  Data[,"idh(Se):units"] <- 0 #this is a dummy column to convince multree that the random effects equation is OK
  Data$O.C <- factor(Data$O.C, labels=c("Stable", "Plastic"))
  return(Data)
}
SetGroups <- function(Dataset){
  Stable <- which(Dataset$O.C == "Stable")
  Plastic <- (1:nrow(Dataset))[-Stable]
  Small <- which(Dataset$SylRepBi == "Smaller")
  SS <- which(Stable %in% Small)
  SL <- Stable[-SS]
  SS <- Stable[SS]
  PS <- Plastic[which(Plastic %in% Small)]
  NewCat <- rep("LargerPlastic", nrow(Dataset))
  NewCat[SL] <- "LargerStable"
  NewCat[SS] <- "SmallerStable"
  NewCat[PS] <- "SmallerPlastic"
  return(NewCat)
}

#Stats Tests and plots
FunnelPlots <- function(Data, Groups="Rep", pos=1){
  colors <- c(rgb(0,0,1,.75), rgb(1,0,0,.75),rgb(0,0,0,.75))
  plot(type='n', Data$Zs,1/Data$Se, xlim=c(-1.5,1.5), ylim=c(min(1/Data$Se),max(1/Data$Se)),
       xlab="Fisher's Transformed Z", ylab = "Inverse Standard Error",
       font=2, font.lab=2, cex.axis=.8, cex.lab=.9)
  abline(v=mean(Data$Zs), lty=3, lwd=2, col="grey40")
  if(Groups=="Rep"){
    log <- log(Data$Syllable.Rep)
    Black <- which(is.na(log))
    minlog <- min(log[-Black])
    maxlog <- max(log[-Black]-minlog)
    colscale <- 1-(log[-Black]-minlog)/maxlog
    points(Data$Zs[-Black], 1/Data$Se[-Black], pch=21, col="black",bg=rgb(1,colscale,colscale,1))
    points(Data$Zs[Black], 1/Data$Se[Black], pch=19)
    scale <- seq(0,1,.05)
    text(-1.35,6.75,"Smaller", adj=0, font=2)
    text(-1.35,9.5,"Larger", adj=0, font=2)
    points(rep(-1.45,length(scale)), 6.75+scale*2.75,
           pch=21, bg=rgb(1,rev(scale),rev(scale),1), col="black")
    #rect(-1.5, 7+scale*2, -1.4,7.5+scale*2, border=NA, col=rgb(1,rev(scale),rev(scale),1))
    #rect(-1.5,7,-1.4,9.5)
    points(-1.45, 6.25,pch=19)
    text(-1.35, 6.25, "Unknown", adj=0, font=2)
    mtext(LETTERS[pos], 2, 1.5, FALSE,max(1/Data$Se), las=2)
    return()
  }else if(Groups=="OC"){
    G1 <- which(Data$O.C == "Stable")
    G2 <- which(Data$O.C == "Plastic")
    LegTex <- c("Song-Stable", "Song-Plastic", "Unknown")
  }else if(Groups=="Metric"){
    G1 <- which(Data$MClass == "Time")
    G2 <- which(Data$MClass == "Offspring")
    LegTex <- c("Time Latency", "NumOffspring", "NumFemale")
  }
  points(Data$Zs[G1], 1/Data$Se[G1], pch=19, col=colors[1])
  points(Data$Zs[G2], 1/Data$Se[G2], pch=19, col=colors[2])
  points(Data$Zs[-c(G1,G2)], 1/Data$Se[-c(G1,G2)], pch=19, col=colors[3])
  par(font=2)
  legend("topleft", LegTex,
         cex=.9,pch=19, pt.bg="black", col=colors)
  par(font=1)
  mtext(LETTERS[pos], 2, 1.5, FALSE,max(1/Data$Se), las=2)
}
AsyncTests <- function(Indicator, Se){
  print(regtest(Indicator, sei=Se))
  print(ranktest(Indicator, sei=Se))
}

#Main MetaAnalysis
SetPrior <- function(Ran, Var=1, Nu=.002, a=1000){
  RanVar <- unlist(strsplit(Ran,"[+]"))
  if("idh(Se):units" %in% RanVar){#vars if measurement error
    RandomN <- length(RanVar)-1
    Add <- list(V = diag(1), fix = 1)
  }else{#vars without measurement error
    RandomN <- length(RanVar)
    Add <- NULL
  }
  Gs <- list(V=matrix(Var/(RandomN+1)), nu=Nu, alpha.mu=rep(.1,1), alpha.V=diag(1)*a)
  prior<-list(R=list(V=matrix(Var/(RandomN+1)),nu=Nu),
              G=rep(list(Gs), RandomN))
  prior$G[[RandomN+1]] <- Add
  return(prior)
}
BayesMeta <- function(Data, Fix="Zs~1", Ran="~1",
                      Vari=1, NU=1, TreeA = NULL, Nitt=200000, Burnin=30000,
                      PR=FALSE){
  #Setup basic prior
  Prior <- SetPrior(Ran,Vari, NU)

    if(!is.null(TreeA)){
    #set up tree
    drop <- which(is.na(match(rownames(TreeA), Data$Species)))
    TreeAinvSubSet <- TreeA[-drop,-drop]
    Meta <- MCMCglmm (as.formula(Fix),random=as.formula(Ran),prior=Prior,
                      ginverse=list(animal=TreeAinvSubSet), #for phylogeny
                      #mev=Data$Se^2,
                      data=Data, nitt=Nitt, burnin=Burnin, verbose=FALSE,
                      pr=PR)
  }else{
    Meta <- MCMCglmm (as.formula(Fix),random=as.formula(Ran),prior=Prior,
                      #mev=Data$Se^2,
                      data=Data, nitt=Nitt, burnin=Burnin, verbose=FALSE,
                      pr=PR)
  }
  return(Meta)
}
BayesMetaMultiTree <- function(Data, MultiTree, Name, Fix="Zs~1", Ran="~animal",
                               Vari=1, NU=1, Nitt=400000, Burnin=30000,
                               PR=FALSE){
  Prior <- SetPrior(Ran,Vari, NU)
  MultiTree <-lapply(MultiTree, function(tree) fixTree(tree, Data$Species))
  class(MultiTree) <- "multiPhylo"
  mulTree.data <- as.mulTree(Data, MultiTree, taxa = "Species",
                             rand.terms = as.formula(Ran))
  mulTree(mulTree.data, as.formula(Fix), parameters = c(Nitt, 10, Burnin),
          priors = Prior, output = Name)
}

fixTree <-function(Tree, Species){
    Tree$tip.label <- gsub("_", " ", Tree$tip.label)
    drop <- which(is.na(match(Tree$tip.label, Species)))
    Tree <- drop.tip(Tree, drop)
    Tree <- multi2di(Tree, .001)
    Tree <- force.ultrametric(Tree)
    Tree$edge.length[which(Tree$edge.length == 0)] <- .0001
    Tree <- force.ultrametric(Tree)
    return()
}

#Basic Graphs
DistributionPlotV2 <- function(Data, lines=TRUE, specsize=.8){
  #pull data, combine, kill dups, sort, separate out closed-ended learners
  Data <- Data[which(is.na(Data$Syllable.Rep)==FALSE),]
  Data <- Data[order(Data$Syllable.Rep),]
  nodups <-Data[which(duplicated(as.character(Data$Species))==FALSE),]
  NumStud <- vector("numeric", length(nodups$Species))
  for(i in seq_along(NumStud)){
    SpecInd <- which(Data$Species == nodups$Species[i])
    NumStud[i] <- length(unique(Data$Study[SpecInd]))
  }
  Col <- rep("black", length(nodups$Species))
  Col[which(nodups$O.C==1 | nodups$O.C== "Stable")] <- "blue"
  Col[which(nodups$O.C==2| nodups$O.C== "Plastic")] <- "red"
  #Close <-nodups[Cind,]
  #stuff for graphing
  ticks <- length(nodups$Species)
  if(lines == TRUE){
    par(mar = c(15,5,.5,1), mgp =c(1.6,.6,0)) 
  }else{
    par(mar = c(15,3,1.5,1), mgp =c(1.6,.6,0))
  }
  
  
  #plot for ln() graph
  Syllable.Rep <- log(nodups$Syllable.Rep)
  #Close$repsize <- log(Close$repsize)
  plot(0, type="n", xlim = c(1,ticks), ylim = c(0, max(Syllable.Rep)),
       las = 2, xaxt = "n", xlab = "", ylab ="ln(Repertoire Size)",  cex = .6,
       cex.axis = .75, font.lab = 2, font.axis=2,
       #Done first so line not crossing axis
       panel.first = {
         if(lines==TRUE){
         rect(which(nodups$Syllable.Rep==17.43)-.5,-1,
              which(nodups$Syllable.Rep==216)-.5,8,
              col = "grey90", border=NA)
         rect(which(nodups$Syllable.Rep==18.15)-.5,-1,
              which(nodups$Syllable.Rep==22.5)-.5,8,
              col = "grey75", border=NA)
           abline(v=c(8.5, 15.5, 21.5), col="black", lty=2)
         }
         })
  
  segments(2:ticks-1, Syllable.Rep[2:ticks-1], 2:ticks, Syllable.Rep[2:ticks],
           lty = 1, lwd = 1.5, col = "black")
  points(1:ticks, Syllable.Rep, pch=19, col = Col, cex = 1+(NumStud-1)*.4)
  axis(1, 1:ticks, labels = paste(nodups$Species, "  (",nodups$Syllable.Rep,")", sep=""), cex.axis = specsize, las=2, font  = 4)
  par(font=2)
  legend("topleft", legend=c("Song-Stable", "Song-Plastic  ", "Unknown"),
         col=c("blue", "red", "Black"), pch = 19, cex=.8)
  par(font=1)
  mtext(side=1, "Species", line=12, font=2)
}
MinMaxDisPlot<- function(DataSet){
  DataSetMAX <- DataSet
  DataSetMAX$Syllable.Rep <- DataSet$Syl.Max
  DistributionPlotV2(DataSetMAX, FALSE, .6)
  mtext("Max Values", font=2)
  mtext("A", side=2, las=2, at=7, line=1.5, font=2)
  
  DataSetMIN <- DataSet
  DataSetMIN$Syllable.Rep <- DataSet$Syl.Min
  DistributionPlotV2(DataSetMIN, FALSE, .6)
  mtext("Min Values", font=2)
  mtext("B", side=2, las=2, at=7, line=1.5, font=2)
}

ThresholTestingPlot <- function(ThreshData, CUT=38){
  Thresh <- sort(unique(ThreshData$Syllable.Rep))
  Thresh <- Thresh[2:(length(Thresh)-1)]
  pvals <- data.frame(Small = numeric(), Large = numeric())
  Means <- pvals
  Cred95 <- data.frame(SmallL = numeric(), SmallU = numeric(),
                       LargeL = numeric(), LargeU = numeric())
  Cred50 <- Cred95
  for(i in seq_along(Thresh)){
    set.seed(49)
    ThreshData[,"SylRepBi"] <- cut(ThreshData$Syllable.Rep, c(0,Thresh[i],2000), labels=c("Smaller", "Larger"), right=FALSE)
    modelThresh <- BayesMeta(Data=ThreshData, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = p.var, TreeA=TreeAinv)
    summary(modelThresh)
    Means[i,] <- c(mean(modelThresh$Sol[,1]), mean(modelThresh$Sol[,2]))
    Cred95[i,] <- c(coda::HPDinterval(modelThresh$Sol[,1])[1:2], coda::HPDinterval(modelThresh$Sol[,2])[1:2])
    Cred50[i,] <- c(coda::HPDinterval(modelThresh$Sol[,1], prob=.5)[1:2],  coda::HPDinterval(modelThresh$Sol[,2], prob=.5)[1:2])
    pvals[i,] <- 2*pmax(0.5/dim(modelThresh$Sol)[1], pmin(colSums(modelThresh$Sol[,, drop = FALSE] > 0)/dim(modelThresh$Sol)[1], 
                                                          1 - colSums(modelThresh$Sol[,, drop = FALSE] > 0)/dim(modelThresh$Sol)[1]))
  }
  #cbind(Thresh,pvals)
  set1 <- 1:(length(Thresh)-1)
  set2 <- 2:length(Thresh)
  set3 <- log(Thresh[set1])
  set4 <- log(Thresh[set2])
  set5 <- 1:length(Thresh)
  set6 <- log(Thresh)
  pdf("TestIntervals.pdf")
  par(mfrow=c(1,2), mar=c(3,3,.5,.5), mgp=c(1.5,.5,0))
  plot(1,type="n", xlim=c(min(log(Thresh)), max(log(Thresh))), ylim=c(.001,1), log='y',
       ylab="pvalue", xlab="ln(Repertoire Size)")
  pMCMCFull <- pMCMC <- 2*pmax(0.5/dim(model8.1$Sol)[1], pmin(colSums(model8.1$Sol[,, drop = FALSE] > 0)/dim(model8.1$Sol)[1], 
                                                              1 - colSums(model8.1$Sol[,, drop = FALSE] > 0)/dim(model8.1$Sol)[1]))
  abline(h=pMCMCFull, lty=2, col="black")
  abline(v=log(CUT))
  segments(set3,pvals[set1,1], set4,pvals[set2,1], lwd=2, col="blue")
  segments(set3,pvals[set1,2], set4,pvals[set2,2], lwd=2, col="red")
  
  plot(1,type="n", xlim=c(min(log(Thresh)), max(log(Thresh))), ylim=c(-1,1.5),
       ylab="Posterior Mean", xlab="ln(Repertoire Size)")
  abline(h=mean(model8.1$Sol[,1]), lty=2, col="black")
  abline(h=0, lwd=2)
  abline(v=log(CUT))
  polygon(x=c(set6, rev(set6)), y=c(Cred95[set5,1],rev(Cred95[set5,2])), col=rgb(0,0,1,.1), border=rgb(0,0,1,.5))
  polygon(x=c(set6, rev(set6)), y=c(Cred95[set5,3],rev(Cred95[set5,4])), col=rgb(1,0,0,.1), border=rgb(1,0,0,.5))
  polygon(x=c(set6, rev(set6)), y=c(Cred50[set5,1],rev(Cred50[set5,2])), col=rgb(0,0,1,.5), border=rgb(0,0,1,.4), density=50)
  polygon(x=c(set6, rev(set6)), y=c(Cred50[set5,3],rev(Cred50[set5,4])), col=rgb(1,0,0,.5), border=rgb(1,0,0,.4), density=50)
  segments(set3, Means[set1,1], set4,Means[set2,1], lwd=2, col="blue")
  segments(set3, Means[set1,2], set4,Means[set2,2], lwd=2, col="red")
}
BayesForestPlot <- function(Data, title="Forest Plot", Group="Rep",
                            mainfont =.7, model1=NULL, model2=NULL){
  #if(is.null(model1) && is.null(model2)){
  #  layout(matrix(c(rep(1,6),2,3,4,5,6,7), nrow=2, ncol=6,byrow=TRUE),
  #         widths = c(.22,.06,.1,.1,.15,.36), heights = c(.05, .95)) 
    layout(matrix(c(rep(1,6),2,3,4,5,6, 7, rep(8,6)), nrow=3, ncol=6,byrow=TRUE),
           widths = c(.22,.06,.1,.1,.15,.36), heights = c(.03, 1, .07)) 

  if(Group == "OC"){
    Data <- Data[order(Data$O.C,decreasing = FALSE),]
    BoxCol <- ifelse(Data$O.C=="Closed", "blue",ifelse(Data$O.C=="Plastic", "red", "black"))
    BottomLab <- c("Song-Stable", "Song-Plastic")
  }else if(Group == "Rep"){
    Data <- Data[order(Data$Syllable.Rep, decreasing = FALSE),]
    BoxCol <- "red"
    BottomLab <- c("Smaller", "Larger")
  } else{
      BoxCol <- "grey50"
  }
  Start <- length(Data$Species)
  PlotSize <- Start+2
  #Title
  par(mar=rep(.01,4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(.3,.3,title, adj=.3, cex=1.2, font=2)
  
  par(mar=c(3,.01,.01,.01))
  
  #species/Study labels
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(1,PlotSize))
  text(x=1,y=PlotSize, "Species (Rep) [citation]", adj=1,font=4)
  text(x=1, y=Start:1,
       paste0(Data$Species, " (", Data$Syllable.Rep, ") ", Data$Study), adj=1, font=3, cex=mainfont)
  #nBirds Labels
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(1,PlotSize))
  text(x=.5,y=PlotSize, "nBirds", adj=.5, font=2)
  text(x=.5, y=Start:1, Data$nBirds, adj=.5, cex=mainfont)
  #Weight Labels
  #plot.new()
  #plot.window(xlim=c(0,1), ylim=c(1,length(PairDate$species.name)+1))
  #text(x=.5,y=length(PairDate$species.name)+1, "Weight", adj=.5, font=2)
  #text(x=.5, y=length(PairDate$species.name):1, paste("idk", "%", sep=""), adj=.5)
  
  #Correlation Labels
  #plot.new()
  #plot.window(xlim=c(0,1), ylim=c(1,PlotSize))
  #text(x=.5,y=PlotSize, "Raw Cor", adj=.5, font=2)
  #text(x=.5, y=Start:1, Data$Raw.Cor, adj=.5, cex=mainfont)
  
  #Correlation Labels
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(1,PlotSize))
  text(x=.5,y=PlotSize, "Fisher's Z", adj=.5, font=2)
  text(x=.5, y=Start:1, round(Data$Zs, digits=3), adj=.5, cex=mainfont)
  
  #Confidence Interval
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(1,PlotSize))
  upper <- round(Data$Zs - qnorm(.05/2, lower.tail = FALSE) * Data$Se, digits=3)
  lower <- round(Data$Zs + qnorm(.05/2, lower.tail = FALSE) * Data$Se, digits=3)
  text(x=.5,y=PlotSize, "[95% CI]", adj=.5, font=2)
  text(x=.5, y=Start:1,
       paste("[",upper,";",lower,"]"), adj=.5, cex=mainfont)
  
  #Forest Plot
  par(mar=c(3,.01,2,.01), mgp=c(1.5,.5,0))
  plot(0,type='n',xlim=c(-1.5,1.5), ylim=c(1,PlotSize),
       bty='n', yaxt='n',
       xlab="Fisher's Z", font.lab=2)
  abline(v=0, lwd=2)
  abline(v=mean(Data$Zs), lty=2, lwd=2, col="grey50")
  bottom <- 1
  Fill <- .5
  rect(xright=Data$Zs-.1,xleft=Data$Zs+.1,
       ybottom=seq(PlotSize+Fill, bottom, length.out = Start),
       ytop=seq(PlotSize+Fill, bottom, length.out = Start)+.5,
       col=BoxCol, border= NA)
  segments(x0=Data$Zs,x1=Data$Zs,
           y0=seq(PlotSize+Fill, bottom, length.out = Start)+.125,
           y1=seq(PlotSize+Fill, bottom, length.out = Start)+.375)
  segments(x0=upper,x1=lower,
           y0=seq(PlotSize+Fill, bottom, length.out = Start)+.25,
           y1=seq(PlotSize+Fill, bottom, length.out = Start)+.25)
  if(Group == "Rep"){
    place <- sapply(c(18.5, 38, 216),
                    function(x) min(which(sort(Data$Syllable.Rep) == x)))
    #print((seq(PlotSize+Fill, 0, length.out = Start))[place-1])
    abline(h=c(71.7,45.9,9), lty=2, col="grey50")
  }
  plot.new()
  
  #if(is.null(model1)==FALSE){
  #  par(mar=c(.01,.01,.01,.01))
  #  plot.new()
  #  plot.window(xlim=c(0,1), ylim=c(0,1))
  #  Printout <- c(mean(model1$Sol[,1]),coda::HPDinterval(model1$Sol),
  #                effectiveSize(model1$Sol[, 1, drop = FALSE]),
  #                2*pmax(0.5/dim(model1$Sol)[1], pmin(colSums(model1$Sol[,1, drop = FALSE] > 0)/dim(model1$Sol)[1], 
  #                                                  1 - colSums(model1$Sol[, 1, drop = FALSE] > 0)/dim(model1$Sol)[1])) )
  #  Printout <- round(Printout, digits=4)
  #  if(Printout[5] == 0){
  #    Printout[5] <- "<0.0001"
  #  }
   # text(seq(.05,.6, length.out=6), .95, c("Group", "post.mean", "l-95% CI", "u-95% CI", "eff.samp",  "pMCMC"), font=2, cex=mainfont)
  #  text(seq(.05,.6, length.out=6), .70, c("Full Population",Printout), cex=mainfont)
  #  if(is.null(model2)==FALSE){
  #    i <- 2
  #    for(i in seq_along(colnames(model2$Sol))){
  #      Printout <- data.frame(mean(model2$Sol[,i]),coda::HPDinterval(model2$Sol[,i]),
  #                    effectiveSize(model2$Sol[, i, drop = FALSE]),
  #                    2*pmax(0.5/dim(model2$Sol)[1], pmin(colSums(model2$Sol[,i, drop = FALSE] > 0)/dim(model2$Sol)[1], 
  #                                                    1 - colSums(model2$Sol[, i, drop = FALSE] > 0)/dim(model2$Sol)[1])) )
  #      Printout <- round(Printout, digits=4)
  #      Printout[4] <- ceiling(Printout[4])
  #      if(Printout[5] == 0){
  #        Printout[5] <- "<0.0001*"
  #      }else if(Printout[5] < .025){
  #        Printout[5] <- paste(Printout[5], "*")
  #      }
  #      text(seq(.05,.6, length.out=6), .7-(.25*i), c(BottomLab[i],Printout), cex=mainfont)
  #    }
  #  }
  #}
}

#Misc Data pulling
ModelDataPuller <- function(ModelList, Group){
  counter <- 1
  DataFrame <- data.frame(character(),character(),
                          character(),character(),
                          stringsAsFactors = FALSE)
  names(DataFrame) <- c("Group", "Post Mean", "95% CredInt", "pMCMC")
  for(i in seq_along(ModelList)){
    for(j in seq_along(Group)){
      credInt <- coda::HPDinterval(ModelList[[i]]$Sol[,j])
      DataFrame[counter,] <- cbind(Group[j],
                                   round(mean(ModelList[[i]]$Sol[,j]), digits=3),
                                   paste0("[", round(credInt[1], digits=3), ";", round(credInt[2], digits=3), "]"),
                                   round(2*pmax(0.5/dim(ModelList[[i]]$Sol)[1], pmin(colSums(ModelList[[i]]$Sol[,j, drop = FALSE] > 0)/dim(ModelList[[i]]$Sol)[1], 
                                                                                     1 - colSums(ModelList[[i]]$Sol[, j, drop = FALSE] > 0)/dim(ModelList[[i]]$Sol)[1])), digits=3))
      counter <- counter+1
    }
  }
  return(DataFrame)
}
BESTPrinter <- function(Dataset, Group, Threshs=NULL, Digits=3, PCrit=2.5){
  ColNames <- c("BEST Mean", "95% CredInt", "%<0")
  DataFrame <- data.frame(character(),character(),character(),
                          stringsAsFactors=FALSE)
  names(DataFrame) <- ColNames
  if(is.null(Threshs)){
    DataFrame[1,] <- BESTSpooler(Dataset, Group, Digits, PCrit)
  }else{
    for(i in seq_along(Threshs)){
      
      Dataset$SylRepBi <-  cut(Dataset$Syllable.Rep, c(0,Threshs[i],2000), labels=c("Smaller", "Larger"), right=FALSE) 
      DataFrame[i,] <- BESTSpooler(Dataset, Group, Digits, PCrit)
    }
    DataFrame[,"Threshold"] <- Threshs
    DataFrame <- DataFrame[,c(4,1:3)]
  }
  return(DataFrame)
}
BESTSpooler <- function(Dataset, Group, Digits=3, PCrit=2.5){
  set.seed=49
  Groups <- levels(Dataset[,Group])
  if(length(Groups) == 0){
    Groups <- unique(Dataset[,Group])
  }
  if(Groups[1] != "Increasing"){
    Groups <- rev(Groups)
  }
  SL <- BESTmcmc(Dataset$Zs[which(Dataset[,Group] == Groups[1])],
                 Dataset$Zs[which(Dataset[,Group] == Groups[2])],
                 numSavedSteps=16998, rnd.seed=49)
  paramSampleVec <- SL$mu1-SL$mu2
  by <- diff(hdi(paramSampleVec))/18
  breaks <- unique(c(seq(from = min(paramSampleVec), 
                         to = max(paramSampleVec), by = by), max(paramSampleVec)))
  histinfo <- hist(paramSampleVec, breaks = breaks, plot = FALSE)
  cvHt <- 0.7 * max(histinfo$density)
  pcgtCompVal <- round(100 * sum(paramSampleVec < 0)/length(paramSampleVec), 
                       1)
  if(pcgtCompVal < PCrit){
    if(pcgtCompVal == 0){
      pcgtCompVal <- "0.1*"
    }else{pcgtCompVal <- paste0(pcgtCompVal,"*")}
  }
  HDI <- hdi(paramSampleVec, .95)
  
  data <- cbind(round(mean(HDI), digits=Digits),
                paste0("[", round(HDI[1], digits=Digits), ";", round(HDI[2], digits=Digits), "]"),
                pcgtCompVal)
  return(data)
}
MetaTable <- function(DataSet, Model, Group, Measure=TRUE, Species="Species",PSort=FALSE, Pcrit = .025, Thresher=NULL, Groups=NULL){
  for(i in 1:length(Model)){
    if(is.null(Groups)){
      reset <- TRUE
      Groups <- levels(DataSet[[i]][,Group])
      if(length(Groups) == 0){
        Groups <- unique(DataSet[[i]][,Group])
      }
      if(Groups[1] == "Stable"){
        Groups <- rev(Groups)
      }
      if(Groups[1] == "Not Increasing"){
        Groups <- rev(Groups)
      }
    }else{reset<-FALSE}
    #Groups <- sort(Groups)
    Base <- ModelDataPuller(list(Model[[i]]), Groups)
    if(Measure){
      Table <- table(as.character(DataSet[[i]][,Group]))
      set <- numeric(length(Groups))
      
      if(DataSet[[i]][1,Group] %in% Groups){
        for(j in seq_along(Groups)){
          set[j] <- Table[which(names(Table) == Groups[j])]
        }
      }else{set <- rep(Table, length(set))}
      Base["#Measure"] <- set
      Base[,"#Measure"] <- as.integer(Base[,"#Measure"])
      Base <- Base[,c(1, ncol(Base), 2:(ncol(Base)-1))]
    }
    if(!is.na(Species)){
      if(DataSet[[i]][,Group][1] %in% Groups){
        Levels <- Groups
      }else{
        Levels <- unique(DataSet[[i]][,Group])
      }
      Base["#Species"] <- sapply(1:length(Levels), function(x) length(unique(
        DataSet[[i]][,Species][
          which(DataSet[[i]][,Group] == Levels[x])])))
      
      Base <- Base[,c(1, ncol(Base), 2:(ncol(Base)-1))]
    }
    if(PSort){
      Base <- Base[order(Base$pMCMC),]
    }
    if(i==1){
      Final <- Base
    }else{
      Final <- rbind(Final, Base)
    }
    if(reset){
      Groups <- NULL
    }
  }
  stars <- which(Final$pMCMC < Pcrit)
  if(length(stars)>0){
    Final$pMCMC <- as.character(Final$pMCMC)
    Final$pMCMC[stars] <- paste0(Final$pMCMC[stars], "*")
  }
  if(!is.null(Thresher)){
    Final[,"Threshold"] <- paste0(c("<", ">="),rep(Thresher, each=2))
    Final <- Final[c(ncol(Final), 1:(ncol(Final)-1))]
    
  }
  return(Final)
}
VarianceTable <- function(ModelList, SE, OutputLatex=TRUE, Caption="", Digits = 2){
  #Make Table
  Table <- data.frame(Fixed=character(),Random=character(),
                      I2=character(), DIC=character(),
                      stringsAsFactors = FALSE)
  count <- 1
  mVar <- (sum(1/SE^2)*(length(SE)-1))/(sum(1/SE^2)^2 - sum((1/SE^2)^2))
  for(i in seq_along(ModelList)){
    VarsVector <- getVars(ModelList[[i]])
    TotalVar <- sum(VarsVector)+mVar
    #get variance data, Get first row of data
    VarNames <- colnames(ModelList[[i]]$VCV)
    VarNames <- gsub("animal", "Phylo", VarNames)
    if(VarsVector[1] <=0){
      FinalFormat <- "Not Well Sampled"
    }else{
      FinalFormat <- paste0(round(VarsVector[1]/TotalVar*100, digits=Digits), "%")
    }
    Table[count,] <- c(deparse(ModelList[[i]]$Fixed$formula),
                       VarNames[1],
                       FinalFormat,
                       round(ModelList[[i]]$DIC, digits=Digits))
    #if more than one random variable (ignoring units and SE
    NumVar <- ncol(ModelList[[i]]$VCV)-1-length(grep("[.]units", VarNames))
    if(NumVar > 1){
      for(j in 2:NumVar){
        if(VarsVector[j] <=0){
          FinalFormat <- "Not Well Sampled"
        }else{
          FinalFormat <- paste0(round(VarsVector[j]/TotalVar*100, digits=Digits), "%")
        }
        Table[count+(j-1),] <- list('""',
                                    VarNames[j],
                                    FinalFormat,
                                    '""') 
      }
    }
    count <- count+NumVar
  }
  if(OutputLatex){
    QuickLatexTable(Table, Caption)
  }else{return(Table)}
}

getVars <- function(model){
  VarData <- numeric(ncol(model$VCV)-2)
  for(k in seq_along(VarData)){
    dense <- density(model$VCV[,k])
    VarData[k] <- dense$x[which.max(dense$y)]
  }
  return(VarData)
}

QuickLatexTable <- function(Table, Caption, Path=file.path(getwd(),"Tables.txt")){
  if(!file.exists(Path)){
    file.create(Path)
  }
  sink(file.path(Path), TRUE)
  print(xtable(Table,caption=Caption), include.rownames=FALSE, caption.placement = "top")
  sink(NULL)
}

#COnvergence plots
ConvergePlot <- function(model, Name="", Alpha.Cex=1, PDF=TRUE){
  if(Name != ""){
    Dash <- "-"
  }else{Dash <- ""}
  
  if(PDF){
    pdf(paste0("Convergence",Dash,Name,".pdf"))
    
  }
  ncols <- ceiling((ncol(model$Sol)+ncol(model$VCV))/2)
  if(ncols > 4){ncols <- 4}
  par(mfrow=c(ncols,4), mgp=c(1.5,.5,0), mar=c(3,3,2,3))
  Count <- 1
  if(ncol(model$Sol)>1 && colnames(model$Sol)=="(Intercept)"){
    Fix <- FALSE    
  }else(Fix <- TRUE)
  for(i in 1:ncol(model$Sol)){
    CheckPlot(model$Sol[,i], colnames(model$Sol)[i], Count, Alpha.cex=Alpha.Cex, Fix)
    Count <- Count+1
    Fix <- TRUE
  }
  for(i in 1:ncol(model$VCV)){
    CheckPlot(model$VCV[,i], colnames(model$VCV)[i], Count, Alpha.cex=Alpha.Cex)
    Count <- Count+1
  }
  if(PDF){
    dev.off()
  }
}
makeChains<- function(model1, Title, Data, Fix, Ran, Vari, TreeA, PR=FALSE){
  model2 <- BayesMeta(Data=Data, Fix=Fix, Ran=Ran, Vari = Vari, TreeA=TreeA, PR=PR)
  model3 <- BayesMeta(Data=Data, Fix=Fix, Ran=Ran, Vari = Vari, TreeA=TreeA, PR=PR)
  getRHat(list(model1, model2, model3), Title)
}
getRHat <- function(ModelList, text="model"){
  sink("RHat.txt", TRUE)
  m <- list(ModelList[[1]]$Sol, ModelList[[2]]$Sol, ModelList[[3]]$Sol)
  print(text)
  print(gelman.diag(do.call(mcmc.list, m)))
  print("")
  gelman.plot(do.call(mcmc.list, m), autoformat=FALSE)
  sink(NULL)
}
CheckPlot <- function(modelVec, title="", pos=1, Alpha.cex=2.5, Fix=TRUE){
  if(length(grep("[.]units",title))>0){
    title <- gsub("[.]units", "", title)
    Se<-TRUE
  }else{Se <- FALSE}
  title <- gsub("[[:punct:]]", "", title)
  if(Fix){
    title <- FixTitles(title)
  }
  plot.default(modelVec, type='l', ylab="", xlab= "Iterations",
               main=paste0("Trace of ", title), cex.main=1)
  if(Se){
    mtext(LETTERS[2*pos-1], 2,1.5, las=2, FALSE, 1.4,cex=Alpha.cex)
  }else{
    mtext(LETTERS[2*pos-1], 2,1.5, las=2, FALSE, max(modelVec), cex=Alpha.cex)
  }
  dense <- density(modelVec)
  plot(dense, cex.main=1,
       main=paste0("Density of ", title), font.main=2)
  segments(modelVec, 0, modelVec, max(dense$y)*.03)
  mtext(LETTERS[2*pos], 2,1.5, las=2, FALSE, max(dense$y),cex=Alpha.cex)    
}
FixTitles <- function(Title){
  if(Title == "animal"){
    Title <- "Phylogenic Effects"
  }else if(Title=="Se"){
    Title <- "Standard Error"
  }else if(Title == "MType"){
    Title <- "Measure Type"
  }else if(Title =="Intercept"){
    Title <- "Population"
  }else if(Title=="Species"){
    Title <- "Species Effect"
  }else if(Title == "units"){
    Title <- "Units"
  }else if(Title=="OCPlastic"){
    Title <- "Plastic"
  }else if(Title =="OCStable"){
    Title = "Stable"
  }else if(Title =="SylRepBiLarger"){
    Title = "Larger"
  }else if(Title =="SylRepBiSmaller"){
    Title = "Smaller"
  }else if(Title == "lnRep"){
    Title = "Ln(Repertoire Size)"
  }else if(Title=="SylRepBiSmallerOCStable"){
    Title = "Smaller Stable"
  }else if(Title=="SylRepBiLargerOCStable"){
    Title = "Larger Stable"
  }else if(Title=="SylRepBiSmallerOCPlastic"){
    Title = "Smaller Plastic"
  }else if(Title=="SylRepBiLargerOCPlastic"){
    Title = "Larger Plastic"
  }
  return(Title)
}

#Extra Analyses
JackThresh <- function(dataset, Vari, Tree, Type="OCJack", Group = "OC"){
  #changes the OC grouping of a study or species
  #aternatively, tests all thresholds
  
  #Set up Data Structures
  if(Type == "Thresh"){
    Thresh <- sort(unique(dataset$Syllable.Rep))
    Thresh <- Thresh[2:(length(Thresh)-1)]
    NumItt <- length(Thresh)
    
    pvals <- data.frame(Small = numeric(NumItt), Large = numeric(NumItt), row.names=Thresh)
    Means <- pvals
    Cred95 <- data.frame(SmallL = numeric(NumItt), SmallU = numeric(NumItt),
                         LargeL = numeric(NumItt), LargeU = numeric(NumItt), row.names=Thresh)
    Cred50 <- Cred95
  }else{
    Species <- as.character(unique(dataset$Species))
    NumItt <- length(Species)
    if(Group == "Rep"){
      pvals <- data.frame(Small = numeric(NumItt), Large = numeric(NumItt), row.names=Species)
      Means <- pvals
      Cred95 <- data.frame(SmallL = numeric(NumItt), SmallU = numeric(NumItt),
                         LargeL = numeric(NumItt), LargeU = numeric(NumItt), row.names=Species)
      Cred50 <- Cred95
      FIXED <- "Zs~SylRepBi-1"
      if(Type=="OCJack"){
        warning("This thing you are doing makes no sense.  Changing the OC status will not affect the repsize results.")
      }
    }else{
      pvals <- data.frame(Closed = numeric(NumItt), Open = numeric(NumItt), row.names=Species)
      Means <- pvals
      Cred95 <- data.frame(ClosedL = numeric(NumItt), ClosedU = numeric(NumItt),
                         OpenL = numeric(NumItt), OpenU = numeric(NumItt), row.names=Species)
      Cred50 <- Cred95
      FIXED <- "Zs~O.C-1"
    }
  }
  #OC Jacknife
  if(Type == "OCJack"){
    Jack <- cbind(replicate(length(Species), dataset$O.C))
    JackKnife <- Replacer(dataset, Jack)
    models <- list()
    for(i in seq_along(Species)){
      set.seed(49)
      dataset$O.C <- JackKnife[,i]
      #This is model 8/9.3 in the main analysis!
      model <- BayesMeta(Data=dataset, Fix=FIXED,
                             Ran="~MType+Species+animal+Study+idh(Se):units",
                             Vari = Vari, TreeA=Tree)
      models[[i]] <- model
      #Means[i,] <- c(mean(model$Sol[,1]), mean(model$Sol[,2]))
      #Cred95[i,] <- c(coda::HPDinterval(model$Sol[,1])[1:2], coda::HPDinterval(model$Sol[,2])[1:2])
      #Cred50[i,] <- c(coda::HPDinterval(model$Sol[,1], prob=.5)[1:2],  coda::HPDinterval(model$Sol[,2], prob=.5)[1:2])
      #pvals[i,] <- 2*pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1], 
      #                                                      1 - colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1]))
    }
    #Species Jackknife
  }else if(Type == "SpJack"){
    HoldBack <- list()
    for(i in seq_along(Species)){
      HoldBack[[i]] <- which(dataset$Species == Species[i])
    }
    models <- list()
    for(i in seq_along(Species)){
      set.seed(49)
      Spec <- dataset[-HoldBack[[i]],]
      #This is model 8/9.2 or 8/9.3 in the main analysis!
      model <- BayesMeta(Data=Spec, Fix=FIXED,
                         Ran="~MType+Species+animal+Study+idh(Se):units",
                         Vari = Vari, TreeA=Tree)
      models[[i]] <- model
      #Means[i,] <- c(mean(model$Sol[,1]), mean(model$Sol[,2]))
     # Cred95[i,] <- c(coda::HPDinterval(model$Sol[,1])[1:2], coda::HPDinterval(model$Sol[,2])[1:2])
      #Cred50[i,] <- c(coda::HPDinterval(model$Sol[,1], prob=.5)[1:2],  coda::HPDinterval(model$Sol[,2], prob=.5)[1:2])
      #pvals[i,] <- 2*pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1], 
      #                                                1 - colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1]))
    }
    #Threshold
  }else{
    Thresh <- sort(unique(dataset$Syllable.Rep))
    Thresh <- Thresh[2:(length(Thresh)-1)]
    ThreshData <- dataset
    models <- list()
    for(i in seq_along(Thresh)){

      set.seed(49)
      #This is model 8/9.2 in the main analysis!
      ThreshData[,"SylRepBi"] <- cut(dataset$Syllable.Rep, c(0,Thresh[i],2000), labels=c("Small", "Large"), right=FALSE)
      model <- BayesMeta(Data=ThreshData, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = Vari, TreeA=TreeAinv)
      makeChains(model, Title=paste("SylRep-", Thresh[i]), Data=ThreshData, Fix="Zs~SylRepBi-1", Ran="~MType+Species+animal+Study+idh(Se):units", Vari = Vari, TreeA=TreeAinv)
      
      #Means[i,] <- c(mean(model$Sol[,1]), mean(model$Sol[,2]))
      #Cred95[i,] <- c(coda::HPDinterval(model$Sol[,1])[1:2], coda::HPDinterval(model$Sol[,2])[1:2])
      #Cred50[i,] <- c(coda::HPDinterval(model$Sol[,1], prob=.5)[1:2],  coda::HPDinterval(model$Sol[,2], prob=.5)[1:2])
      #pvals[i,] <- 2*pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1], 
      #                                                    1 - colSums(model$Sol[,, drop = FALSE] > 0)/dim(model$Sol)[1]) )
      models[[i]] <- model
    }
  }
  #Output <- list(Means, Cred95, Cred50, pvals)
  #names(Output) <- c("Posterior.Mean", "Credibility95%", "Credibility50%", "pVals")
  #Output <- list(Output, models)
  return(models)
}
CredibilityIntPlots <- function(JakThr, model, Type){
  pvals <- JakThr$pVals
  Means <- JakThr$Posterior.Mean
  Cred95 <- JakThr$`Credibility95%`
  Cred50 <- JakThr$`Credibility50%`
  pMCMCFull<- 2*pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,1:ncol(model$Sol), drop = FALSE] > 0)/dim(model$Sol)[1], 
                     1 - colSums(model$Sol[,1:ncol(model$Sol), drop = FALSE] > 0)/dim(model$Sol)[1]))
  
  
  #breaks for segment and polygon
  set1 <- 1:(length(rownames(pvals))-1)
  set2 <- set1+1
  set5 <- 1:length(rownames(pvals))
  if(length(model$Sol[1,]) == 1){
    PostMean <- mean(model$Sol)
  }else{
    PostMean <- colMeans(model$Sol)
  }
  #X-axis points 
  if(Type == "Thresh"){
    xaxt <- as.numeric(rownames(pvals)) 
    set3 <- log(xaxt[set1])
    set4 <- log(xaxt[set2])
    set6 <- log(xaxt)
    xlabel <- "Syllable Repertoire Size"
    bmar <- 3
  }else{
    xaxt <- rownames(pvals) 
    set3 <- set1
    set4 <- set2
    set6 <- set5
    if(Type=="SPJack"){
      xlabel <- "Species Removed"
    }else{
      xlabel <- "Species Switched" 
    }
    bmar <- 10
  }
  
  
    
    #pdf("TestIntervals.pdf")
    par(mfrow=c(1,2), mar=c(bmar,4,.5,.5), mgp=c(1.5,.5,0))
    plot(1,type="n", xlim=c(min(set6), max(set6)), ylim=c(.001,1), log='y',
         xaxt="n", ylab="pMCMC", xlab="")
    axis(1, at=set6, labels = xaxt, las=2, cex.axis=.7)
    mtext(xlabel, side=1, line=bmar-1)
    abline(h=pMCMCFull, lty=2, col=c("blue","red"))
    abline(h=.025, lwd=2)
    segments(set3,pvals[set1,1], set4,pvals[set2,1], lwd=2, col="blue")
    segments(set3,pvals[set1,2], set4,pvals[set2,2], lwd=2, col="red")
    
    plot(1,type="n", xlim=c(min(set6), max(set6)), ylim=c(-1,1.5),
         xaxt="n", ylab = "Posterior Mean", xlab="")
    axis(1, at=set6, labels = xaxt, las=2, cex.axis=.7)
    mtext(xlabel, side=1, line=bmar-1)
    abline(h=PostMean, lty=2, col=c("blue","red"))
    abline(h=0, lwd=2)
    polygon(x=c(set6, rev(set6)), y=c(Cred95[set5,1],rev(Cred95[set5,2])), col=rgb(0,0,1,.1), border=NA)
    polygon(x=c(set6, rev(set6)), y=c(Cred95[set5,3],rev(Cred95[set5,4])), col=rgb(1,0,0,.1), border=NA)
    polygon(x=c(set6, rev(set6)), y=c(Cred50[set5,1],rev(Cred50[set5,2])), col=rgb(0,0,1,.5), border=NA, density=50)
    polygon(x=c(set6, rev(set6)), y=c(Cred50[set5,3],rev(Cred50[set5,4])), col=rgb(1,0,0,.5), border=NA, density=50)
    segments(set3, Means[set1,1], set4,Means[set2,1], lwd=2, col="blue")
    segments(set3, Means[set1,2], set4,Means[set2,2], lwd=2, col="red")
}
Replacer <- function(dataset, Jack){
  ind <- which(duplicated(dataset$Species) == FALSE)
  for(i in seq_along(ind)){
    ifelse(ind[i] == "Stable", Switch <- "Plastic", Switch <- "Stable")
    ifelse(i != length(ind), Jack[ind[i]:(ind[i+1]-1),i] <- Switch, Jack[ind[i]:length(Jack[,i]),i] <- Switch) 
  }
  return(Jack)
}

#Arbitrary Model
GenerateArbitraryDataSets <- function(Dataset,Nu=500){
  set.seed(49)
  ArbDataSets <- list()
  Spec <- as.character(unique(Dataset$Species))
  Arbitrary <- Dataset
  for(i in 1:Nu){
    Class <- ifelse(runif(length(Spec)) < .5, "Smaller", "Larger")
    for(j in 1:length(Spec)){
      Arbitrary$SylRepBi[which(Arbitrary$Species == Spec[j])] <- Class[j]
    }
    ArbDataSets[[i]] <- Arbitrary
  }
  return(ArbDataSets)
}
ArbitraryPlots<- function(ArbitraryDataSets, ArbitraryModels, ArbitraryBESTResults, ThresholdTest, Thresholds, Nitt=500){
  #set up original Data
  OriginalMetaData <- as.numeric(ModelDataPuller(ThresholdTest, c("Smaller", "Larger"))$pMCMC)
  OriginalBESTData <- BESTPrinter(SylDataSet, "SylRepBi", Thresholds[2:length(Thresholds)])$`%<0`
  OriginalBESTData <- as.numeric(gsub("[*]", "",OriginalBESTData))
  
  #Set Up arb data
  pMCMC <- MetaTable(ArbitraryDataSets, ArbitraryModels, "SylRepBi")$pMCMC
  pMCMC <- as.numeric(gsub("[*]", "",pMCMC))
  BESTPercentage <- ArbitraryBESTResults$`%<0`
  BESTPercentage <- as.numeric(gsub("[*]", "",BESTPercentage))
  
  #significant models
  SignificantpMCMC <- which(pMCMC<.025)
  SignificantMeta <- ceiling(SignificantpMCMC/2)
  SignificantPercentage <- which(BESTPercentage < 2.5 | BESTPercentage > 97.5)
  SignificantModels <- SignificantMeta[which(SignificantMeta %in% SignificantPercentage)]
  
  par(mar=c(10,3.5,3.8,1),mfrow=c(1,2),mgp=c(1.5,.5,0))
  ArbSub1(pMCMC,OriginalMetaData, BESTPercentage, OriginalBESTData, SignificantModels, Nitt)
  ArdSub2(ArbitraryModels, ArbitraryDataSets, SignificantpMCMC, SignificantMeta)
}
ArbSub1<- function(pMCMCNew, pMCMCOld, NewBEST, OldBEST, SignificantModels, Nitt=500){
  #get the Best values to plot by (must be less than 50%)
  Compare <- cbind(as.numeric(pMCMCNew[seq(1,length(pMCMCNew), by=2)]),
                   as.numeric(pMCMCNew[seq(2,length(pMCMCNew), by=2)]))
  pMCMCXvals <- apply(Compare,1,min)
  Flipper <- which(NewBEST > 50)
  PercentYvals <- NewBEST
  PercentYvals[Flipper] <- ((NewBEST[Flipper])-100)*-1
  
  plot(pMCMCXvals, PercentYvals,
       pch=21, bg=rgb(.4,.4,.4,.5),
       ylim=c(0,50), xlim=c(.001, .1),
       xlab="pMCMC", ylab="BEST %<0", font.lab=2,
       cex=.6, col=rgb(0,0,0,.8),
       panel.first={
         abline(v=c(.067, .013), lwd=2, lty=c(2,3), col=c("black", "grey65"))
         rect(-.2,-.2,.025,2.5, col=rgb(.5,.5,1,.8), border=NA)
       })
  #from the original model
  pMCMCXvals2 <- pMCMCOld[seq(2,length(pMCMCOld), by=2)]
  points(pMCMCXvals2[-1], OldBEST,
         pch=24, bg=rgb(1,0,0,.5), cex=1)
  mtext("A", side=2, las=2, at=52.5, line=1.5, font=2)
  
  vales <- sort(as.numeric(pMCMCNew))
  print(paste0(100*length(SignificantModels)/Nitt,"%"))
  print(vales[floor(500*.025)])
}
ArdSub2 <-function(ArbitraryModels, ArbitraryDataSets, SignificantpMCMC, SignificantMeta){
  Models <- lapply(SignificantMeta, function(x) ArbitraryModels[[x]])
  #MetaTable(ArbitraryDataSets[SignificantMeta], Models, "SylRepBi")
  SigIndex <- SignificantMeta
  LSig <- which(SignificantpMCMC%%2==0)
  x <- 1
  Smaller <- sapply(1:length(Models), function(x) unique(ArbitraryDataSets[[SigIndex[x]]]$Species[
    which(ArbitraryDataSets[[SigIndex[x]]]$SylRepBi=="Smaller")]))
  Larger <- sapply(1:length(Models), function(x) unique(ArbitraryDataSets[[SigIndex[x]]]$Species[
    which(ArbitraryDataSets[[SigIndex[x]]]$SylRepBi=="Larger")]))
  #Large <- sapply(1:length(Models), function(x) unique(Arbitrary$Species[which(Models[[x]]$X[,2]==1)]))
  Species <- unique(as.character(ArbitraryDataSets[[1]]$Species))
  SpecInd <- which(!duplicated(ArbitraryDataSets[[1]]$Species))
  Syls <- ArbitraryDataSets[[1]]$Syllable.Rep[SpecInd]
  SpeciesOrder <- Species[order(Syls)]
  
  TrueLarge <- as.vector(unlist(c(Smaller[-LSig], Larger[LSig])))
  SplitLarge <- sapply(1:length(Species), function(x) length(which(TrueLarge==SpeciesOrder[x])))
  
  
  par(mar=c(11,3.5,0,1),mgp=c(1.5,.5,0))
  hist(rep(1:length(Species),SplitLarge), breaks=0:25, col=rgb(.5,.5,1,.8),
       ylim=c(0, length(SignificantMeta)), xaxt="n", xlab="", main="", yaxt="n",
       ylab="Percentage of Significant Groups With a Species        ",
       font.lab=2, cex.lab=.9)
  
  
  nChoice <- length(SignificantMeta)
  probs <- dbinom(0:nChoice,nChoice,.5)
  Halfway <- cumsum(probs[(ceiling(nChoice/2)+1):nChoice])
  PlusMinus<-max(which(Halfway<.95/2))
  Fifty <- nChoice/2
  rect(0,Fifty-PlusMinus,length(Species),Fifty+PlusMinus, col=rgb(.4,.4,.4,.4), border=NA)
  axis(1, seq(.5,length(Species),by=1), paste0(SpeciesOrder, "(",
                                               sort(Syls),
                                            ")"), las=2, pos=0, cex.axis=.6, font.axis=3)
  axis(2, seq(0,length(SignificantMeta),by=length(SignificantMeta)/10), seq(0,100,by=10), las=2, pos=0)
  mtext("B", side=2, las=2, at=44.5, line=1.5, font=2)
  segments(0,length(SignificantMeta),25,length(SignificantMeta), lwd=1, lty=1)
  segments(25,0,25,length(SignificantMeta), lwd=1, lty=1)
  
  segments(0,length(SignificantMeta)/2,25,length(SignificantMeta)/2, lwd=2, lty=2)
}

