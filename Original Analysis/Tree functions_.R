########
#Coded by Cristina Robinson
#Last Modified 5-1-2018
#Written in R-Studio Version 1.1.442
#R Version 3.4.4
#metafor v2.0-0   Matrix v1.2-12   meta v4.9-1
#phytools v0.6-44
########

require("maps")
require("ape")
require("phytools")

if (file.exists("Output")){
  setwd(file.path(dir, "Output"))
} else {
  dir.create(file.path(dir, "Output"))
  setwd(file.path(dir, "Output"))
  
}

################
#Functions for the rainbow trees
################
LoadPrettyTree <- function(filename){
#what it says: runs the fuctions necessary to create a nice consensus tree
  RawTree=read.nexus(filename)
  ConsensusTree <- consensus.edges(RawTree)
  ConsensusTree <- PrepareConsensus(ConsensusTree)
  ConsensusTree <- PrettyConsensus(ConsensusTree)
  return(ConsensusTree)
}
LoadBirdData <- function(filename, Tree, Pair = FALSE){
#what it says: runs the fuctions necessary to load and process the species data
  birbs <- Loadbirds(filename, Tree)
  birbs <- Dupcleaner(birbs)
  Dataset <- Datasetter(birbs, Pair)
  return(Dataset)
  
}
Loadbirds <- function(filename, concensustree, drop = FALSE, log = FALSE){
#reading in the info  
  birdorder <- c(concensustree$tip.label)
  allinfo <- read.table(filename, sep=",", header=TRUE)
  
##these next two ifs are NOT used in this project!
#reduce large syllable counts to 100
  if(drop == TRUE){
    over100 <- which(allinfo$repsize> 100)
    allinfo$repsize[over100] <- 100 
  }
  if(log == TRUE){
    allinfo$repsize <- log(allinfo$repsize)
  }
  
#reorder to match the concensus tree
  allinfo$species.name <- gsub(" ", "_", allinfo$species.name)
  reorder <- 1:length(allinfo$species.name)
  k <- 1
  i <- 4
  for(i in 1:length(birdorder)){
    ind <- which(allinfo$species.name == birdorder[i])
    Len <- length(ind)
    reorder[k:(k+Len-1)] <- ind
    k <- k+(Len)
  }
  allinforeorder <- allinfo[reorder,]
  
#Convert 2s to 0s
  wayward2 <- which(allinforeorder$o.c == 2)
  allinforeorder$o.c[wayward2] <- 0
  
#not used in this project
#drop the outgroup
  if(drop == TRUE){
    allinforeorder <- allinforeorder[2:length(birdorder),]
  }
  
#Return the final data  
  return(allinforeorder)
}
Datasetter <- function(dataset, Pair = FALSE){
#pull data, binarize the rep size
  repertoiresize <- dataset$repsize
  names(repertoiresize) <- dataset$species.name
  small <- which(repertoiresize < 38)
  repertoiresize[small] <- 0
  repertoiresize[-small] <- 1
  
#name the OC data  
  OpenClose <- dataset$o.c
  OC <- factor(OpenClose, labels = c("open", "close"))
  names(OpenClose) <- dataset$tip.label
  
#name the migratory data  
  Migir <- dataset$migrate
#combines the mixed and migratory birds
  Migir[which(Migir == 1)] <- 2
  names(Migir) <- names(OpenClose)
  
#name the correlations  
  Core <- dataset$cor
  names(Core) <- names(OpenClose)
  if(Pair == TRUE){
    Core <- Core*-1
  }
  
#make a new dataframe that only has the info needed for this analysis
  VarList <- cbind(repertoiresize, OpenClose, Migir, Core)
  return(VarList)
}
Dupcleaner <- function(dataset){
#find duplicated species
  dups <- which(duplicated(dataset$species.name) == TRUE)
  dupspec <- unique(dataset$species.name[dups])
  
#for each duplicated species, find the average correlation
#replace the cor on first entry for that species with the average
  for(i in 1:length(dupspec)){
    ind <- which(dataset$species.name == dupspec[i])
    newcor <- mean(dataset$cor[ind])
    newn <- sum(dataset$nbirds[ind])
    dataset$cor[ind[1]] <- newcor
    dataset$nbirds[ind[1]] <- newn
  }
  
#drop the duplicates
  dataset <- dataset[-dups,]
  return(dataset)
}
PrettyConsensus <- function(conse, sort = FALSE){
  #there is a minor issue with the consensus tree code in that it does not
  #ensure that the full length of each path is the same (so that all of tips end at the same place).
  
  #This code first recreates the node paths from root to tip for each species as a vector of nodes (the loops).
  #Once this exist, it calculates the length of each path and finds the max length.
  #It then subtracts each path length from the max length to get a value of how much
  #length needs to be added to have the tips line up.
  #That value is added to the last edge in the path (the one connecting to the tip).
  
  #Based on my testing, this does not affect future calculations and makes the trees look better.
  #If you are unsure, comment out the "PrettyConsensus" line of LoadPrettyTree;
  #the rest of the code will run just fine. :)
  i <- 0
  numspec <- length(conse[[2]])
  speciespath <- matrix(nrow = numspec, ncol = 3)
  #find new nodes counting forward from i
  for (j in 1:numspec){
    i <- i+1
    startedge <- i
    while(conse[[1]][i,2] > numspec){
      #conse[[1]][i,2] > numspec){
      #Counts up sequentially counting nodes
      while(conse[[1]][i,2] == conse[[1]][i+1,1] && i+1 <= numspec){
        i <- i+1
      }
      #Counts up nonsequentually counting internodes
      if(conse[[1]][i,2] > numspec){
        i <- i+1
      }
    }
    endedge <- i
    nodepath <- startedge:endedge
    
    #if i did not get to the root (i=1) count backwards from i/k to find path
    if(conse[[1]][startedge,1] != numspec + 1){
      k <- startedge
      while(k != 1  && conse[[1]][k,1] != numspec + 1){
        #finds where to jump to in order to go back
        linker <- which(conse[[1]][,2] == conse[[1]][k,1])
        k <- linker
        #count down until reaches a nonsequential number greater than or equal to 1(root)
        
        while(conse[[1]][k,1] == conse[[1]][k-1,2] && k > 1){
          k <- k-1
        }
        nodepath <- c(nodepath, k:linker)
      }
    }
    #get the full path length based on indicies above, get tip index, load into database
    pathleng <- sum(conse[[4]][nodepath])
    finaledge <- conse[[4]][max(nodepath)]
    speciespath[j,]<- cbind(endedge, pathleng, finaledge)
    
  }
  tipstop <- max(speciespath[,2])
  increments <- (speciespath[,2]-tipstop)*-1
  conse[[4]][speciespath[,1]] <- conse[[4]][speciespath[,1]] + increments
  if(sort == TRUE){
    specorder <- conse[[1]][speciespath[,1],2]
    Speciescorrect <- conse$tip.label[specorder]
    return(Speciescorrect)
  } else {return(conse)}
}
PrepareConsensus <-function(conse){
#phytools becomes irate about edge lengths of 0
  conse$edge.length[conse$edge.length==0] = 0.0000001
#the tree needs to be unrooted for phytools.
#We used sayornis phoebe (suboscine) as the root and remove it
#to create our unrooted tree
  conse=drop.tip(conse, 1, root.edge = 0)
  return(conse)
}


DoubleRainbowplot <- function(conse, repsize, OCdata, Migrate, Core,
                              title="RainbowPlot", ANOVA = FALSE, output = TRUE, Fig = ""){
#set up things I need for plotting and sorting by species
  Specieslab <- PrettyConsensus(conse, sort = TRUE)
  Migrate <- Migrate[Specieslab]
  OCdata <- OCdata[Specieslab]
  repsize <- repsize[Specieslab]
  Specieslab <- gsub("_", " ", Specieslab)
  co <- c("white", "black")
  miglab <- c("Sedentary", "Mixed", "Migratory")
  Migrate[which(Migrate == 2)] <- 1

#Get the correlation tree
  CD<-contMap(conse,x=Core,type="fan",fsize=0.05,lwd=2,plot=FALSE)
  names(CD$cols)=rev(names(CD$cols))
  
  if(output == TRUE){
    pdf(paste(title, "Rainbowplot.pdf",sep=" "))
  }

#Start plotting
    layout(matrix(c( rep(1,3), 2:4, rep(5,3)),nrow = 3, byrow=TRUE),
           widths=c(0.42,0.315,0.315), heights=c(.1,.3,.1))
  
#title and column names
  par(mar=c(.1,.1,.1,.1))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(.02,.33, Fig, font = 2, cex = 2)
  n <- c(expression(Correlation~between~Individual~Repertoire),
         bquote(Size~and~bold(.(title))))
  text(.17, c(.19, .07), n, cex = 1.2)
  text(.5+.3, 0.16, "Species", cex = 1.2)

#Legend
  a <- .7
  b <- .23
  text(.438, .85-a, "Large\nSmall", cex = 1, adj = 1)
  text(.560, .85-a, "Closed\nOpen", cex = 1, adj = 1)
  text(.695, .85-a, "Migratory\nSedentary", cex = 1, adj = 1)
  points(rep(.678-b,2), c(.91-a,.79-a), pch=c(24,2), cex = 1.5, bg = 'black')
  points(rep(.8-b,2), c(.91-a,.79-a), pch=c(19,1), cex = 1.5)
  points(rep(.935-b,2), c(.91-a,.79-a), pch=c(15,0), cex = 1.5)

#Plot the correlation data
  par(mar=c(.1,.1,.1,.1))
  CD<-setMap(CD,colors=c("red","white","blue"))
  plot(CD,lwd=3, ftype = "off", legend=FALSE, direction = "rightwards")

#Add in the circles, squares and triangles  
  plot.new()
  plot.window(xlim=c(-0.11,0.05),ylim=c(1, length(conse$tip.label)))
  par(cex=.6)
  points(rep(-.09,length(Specieslab)), 1:length(Specieslab), pch = ifelse(repsize==0, 2, 24), cex = 1.5, bg = 'black')
  segments(-.06,1,-.06,length(conse$tip.label))
  points(rep(-.03,length(Specieslab)), 1:length(Specieslab), pch = ifelse(OCdata==0, 1, 19), cex = 1.5)
  segments(0,1,0,length(conse$tip.label))
  points(rep(.03,length(Specieslab)), 1:length(Specieslab), pch = ifelse(Migrate==0, 0, 15), cex = 1.5)

#Shove the species labels in there
  plot.new()
  plot.window(xlim=c(-0.1,0.1),ylim=c(1, length(conse$tip.label)))
  par(cex=.75)
  text(-.1, 1:length(Specieslab),Specieslab, font = 3, adj = 0)
  
#add figurelegend
  plot.new()
  plot.window(xlim=c(-0.1,0.1),ylim=c(0,1))
  add.color.bar(.03,rev(CD$cols),title="",
                lims=CD$lims,digits=2, prompt=FALSE,x=-.104,
                y=.905,lwd=3,fsize=1, ftype = "b",subtitle="")
  text(-.1, .835, "Correlation", font = 2, adj=0)
  


#tag on the ANOVA data at the end
  if(ANOVA == TRUE){
    N <- 3
#stuff for the text print out
    categor <- c("Repertoire Size", "Learning Style", "Migratory Style")
    Mean1 <- c(mean(Core[which(repsize == 0)]),
             mean(Core[which(OCdata == 0)]),
             mean(Core[which(Migrate == 0)]))
    Mean2 <-c(mean(Core[which(repsize == 1)]),
            mean(Core[which(OCdata == 1)]),
            mean(Core[which(Migrate == 1)]))
    Mean1 <- round(Mean1, digits = 3)
    Mean2 <- round(Mean2, digits = 3)
    names(Mean1) <- c(" Small: ", " Closed: ", " Sedentary: ")
    names(Mean2) <- c(" Large: ", " Open: ", " Migratory: ")
  
#Repsize ANOVA
    print("ANOVA 1")
    RepAnov <- phylANOVA(conse, repsize, Core, nsim = 5000)
#OC ANOVA
    print("ANOVA 2")
    OCAnov <- phylANOVA(conse, OCdata, Core, nsim = 5000)
  
#Migratory ANOVA
    if(length(which(Migrate == 0)) != 0){
      print("ANOVA 3")
      MigAnov <- phylANOVA(conse, Migrate, Core, nsim = 5000)
      Pvals <- c(RepAnov$Pf, OCAnov$Pf, MigAnov$Pf)
      Fval <- c(RepAnov$F, OCAnov$F, MigAnov$F)
    } else {
      N <- N-1
      categor <- categor[1:2]
      Pvals <- c(RepAnov$Pf, OCAnov$Pf)
      Mean1 <- Mean1[1:2]
      Mean2 <- Mean2[1:2]
      Fval <- c(RepAnov$F, OCAnov$F)
    }
    Fval <- round(Fval, digits = 3)  
#Bonferroni correction
    alpha <- .05/N
    alpha <- round(alpha, digits = 3)
  
#Plot the ANOVA text
    par(cex=.6)
    text(-.1,.4,
         paste("Phylogenically-Controlled ANOVA Results:", "\n",
               "Bonferroni-Corrected Alpha = ", alpha, "\n",
               "\n",
               paste(categor, " F-Value = ", Fval,  " p-value = ", Pvals, "  Means:", names(Mean1), Mean1,
                     names(Mean2), Mean2, sep = "", collapse = "\n"), sep = ""), adj = 0)
  }
#Finished!
  if(output == TRUE){
    dev.off()
  }
}

#These funs are used to create the consensus tree if you are using an older version
#of Phytools!
#funs from Liam Revel on Github (liamRevel) and his blog (http://blog.phytools.org/)
#See: liamrevell/phytools/R/utilities.R
#And: liamrevell/phytools/R/consensus.edges.R
if(OLD == TRUE){
  matchLabels<-function(tr1,tr2){
    foo<-function(x,y) if(length(obj<-which(y==x))>0) obj else NA
    M<-cbind(1:Ntip(tr1),sapply(tr1$tip.label,foo,y=tr2$tip.label))
    colnames(M)<-c("tr1","tr2")
    M
  }
  consensus.edges<-function(trees,method=c("mean.edge","least.squares"),...){
    if(hasArg(consensus.tree)) tree<-list(...)$consensus.tree
    else tree<-consensus(trees,p=0.5)
    if(hasArg(if.absent)) if.absent<-list(...)$if.absent
    else if.absent<-"zero"
      N<-length(trees)
    if(method[1]=="mean.edge"){
      M<-lapply(trees,function(x,y) rbind(matchLabels(y,x),matchNodes(y,x)),y=tree)
      nodes<-M[[1]][,1]
      edge.length<-vector(mode="numeric",length=length(nodes))
    for(i in 2:length(nodes)){
      ii<-which(tree$edge[,2]==nodes[i])
      n.absent<-0
      for(j in 1:N){
        edge.length[ii]<-edge.length[ii] +
          if(!is.na(M[[j]][i,2])) trees[[j]]$edge.length[which(trees[[j]]$edge[,2]==M[[j]][i,2])]/N
        else 0
        if(is.na(M[[j]][i,2])) n.absent<-n.absent+1
      }
      if(if.absent=="ignore") edge.length[ii]<-edge.length[ii]*N/(N-n.absent)
    }
    tree$edge.length<-edge.length
  } else if(method[1]=="least.squares"){
    D<-Reduce('+',lapply(trees,function(x,t) cophenetic(x)[t,t],t=tree$tip.label))/N
    tree<-nnls.tree(D,tree=tree,rooted=is.rooted(tree))
  }
  tree
  }
}

