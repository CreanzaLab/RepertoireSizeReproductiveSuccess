########
#Coded by Cristina Robinson
#Last Modified 5-1-2018
#Written in R-Studio Version 1.1.442
#R Version 3.4.4
#metafor v2.0-0   Matrix v1.2-12   meta v4.9-1 
########



require("meta")
require("Matrix")
require("metafor")

if (file.exists("Output")){
  setwd(file.path(dir, "Output"))
} else {
  dir.create(file.path(dir, "Output"))
  setwd(file.path(dir, "Output"))
  
}

################
#Functions for the main meta-analysis
################
QuickMeta <- function(dataset, OC = FALSE, Reper = NULL, Popu = FALSE, Plot = FALSE,
                      PrinPDF = FALSE, ReturnALL = FALSE, ReturnP = FALSE,
                      pred = FALSE, Diag = FALSE){
#code to generate metaanalysis results and a forest plot
#Pull data and sort by rep size s to l
  variablename <- deparse(substitute(dataset))
  dataset <- dataset[order(dataset$repsize),]
  
#Pull vectors of data for convenience
  nbirds <- dataset$nbirds
  studyR <- dataset$cor
  label <- dataset$species.name
  
#set up the group to split on; names the groups and sets the box colors
#Final elseif runs the full dataset with no grouping (population)
  if(is.null(Reper) == FALSE){
    Group <- rep("Smaller Syl Rep", length(dataset$repsize))
    Group[which(dataset$repsize >= Reper)] <- "Larger Syl Rep"
    diagrep <- TRUE
    Repcol <- c(rep("cyan", length(which(Group == "Smaller Syl Rep"))),
                rep("deeppink", length(which(Group == "Larger Syl Rep"))))
  }else if(OC == TRUE){
    Group <- rep("Song-Stable", length(dataset$o.c))
    Group[which(dataset$o.c == 2)] <- "Song-Plastic"
    diagrep <- FALSE
    Repcol <- dataset$o.c
    Repcol[which(Repcol == 1)] <- "blue"
    Repcol[which(Repcol == 2)] <- "red"
  }else if(Popu == TRUE){
    Group <- rep("Population", length(dataset$nbirds))
    diagrep <- NULL
    Repcol <- rep("black", length(Group))
  }

#Labeling the species for the meta-analysis
  subs <- levels(factor(Group))
  if(subs[1] != paste0(factor(Group)[1])){
   subs <- subs[2:1]
  }
  subs <- c(paste("Pooled (", subs, ")", sep = ""))
#run metaanalysis
  analysis <- metacor(studyR, nbirds, label, comb.random=TRUE, comb.fixed=FALSE,
                      byvar = Group, method.tau = "REML", prediction = pred)

#print plot or pdf + generate name
  if(Plot == TRUE){
    forest(analysis, plotwidth=grid::unit(1.5, "cm"), print.I2 = FALSE,
           print.tau2 = FALSE, print.pval.Q = FALSE, col.by = "grey45",
           test.subgroup.random = TRUE, just.addcols = "center", squaresize=.8,
           col.square = Repcol, fontsize = 7, colgap=grid::unit(1.5, "mm"), colgap.forest = "5mm",
           fs.hetstat=7, layout = "revman", addrow = FALSE, spacing = .55,
           xlim=c(-1,1), col.square.lines=Repcol,  digits.pval.Q=4, ff.study.labels = "italic",
           xlab = "COR", text.random = "Full Population",
           text.random.w = subs, print.subgroup.labels = TRUE,
           label.test.subgroup.random = "Subpopulation Differences: ", smlab = "",
           leftlabs = c("Species [ref]", "nBirds", "Weight", "COR [95% CI]"), col.inside="black",
           studlab = paste(dataset$species.name, dataset$citation))
  }else if(PrinPDF == TRUE){
    filename <- paste(variablename, paste(levels(factor(Group)), collapse = ""), Reper, ".pdf", sep="")
    filename <- gsub(" ", "", filename, fixed = TRUE)
    pdf(paste(filename))
    forest(analysis, plotwidth=grid::unit(1.5, "cm"), print.I2 = FALSE,
           print.tau2 = FALSE, print.pval.Q = FALSE, col.by = "grey45",
           test.subgroup.random = TRUE, just.addcols = "center", squaresize=.8,
           col.square = Repcol, fontsize = 7, colgap=grid::unit(1.5, "mm"), colgap.forest = "5mm",
           fs.hetstat=7, layout = "revman", addrow = FALSE, spacing = .55,
           xlim=c(-1,1), col.square.lines=Repcol,  digits.pval.Q=4, ff.study.labels = "italic",
           xlab = "COR", text.random = "Full Population",
           text.random.w = subs, print.subgroup.labels = TRUE,
           label.test.subgroup.random = "Subpopulation Differences: ", smlab = "",
           leftlabs = c("Species", "nBirds", "Weight", "COR [95% CI]"), col.inside="black",
           studlab = paste(dataset$species.name, dataset$citation))
    dev.off()
  }
  
#return metaanlysis data (all, the p-values or I^2 and Tau^2)
  if(ReturnALL == TRUE){
    return(analysis)
  }else if(ReturnP == TRUE){
    p <- Datapuller(analysis)
    return(p)
  }else if(Diag == TRUE){
    ITau <- Datapuller(analysis, type= "diag", Rep = diagrep)
    return(ITau)
  }
}
pvalueplot <- function(dataset, smallestn = 2, replog = TRUE, log = FALSE,
                       INSET = FALSE, pos = 1, short = FALSE, titles = NULL){
#sort data, pull rep, truncate to relevant numbers, make dummy vector
  if(is.null(titles)== TRUE){
  titles <- deparse((substitute(dataset)))
  }
  dataset <- dataset[order(dataset$repsize),]
  fullrep <- dataset$repsize
  nodups <- RepCutter(dataset$repsize, smallestn)
  boundryrep <- fullrep[which(fullrep >= min(nodups) & fullrep <= max(nodups))]
  pstack <- 1:length(nodups)
  
#run quickmeta, spliting on all reps (nodups) to get p-vals
  for(i in 1:length(nodups)){
    pvalue <- QuickMeta(dataset, Reper = nodups[i], ReturnP = TRUE)
    pstack[i] <- pvalue
  }
  
#get the N for small rep at each threshhold value
  ntext <- Ncounter(fullrep, boundryrep)
  if(replog == TRUE){
    nodups <- log(nodups)
  }
#this sets spacing for the plots
  if(log == TRUE){
    texty <- pstack*1.7
    top <- 1.7
    logger <- "y"
    axer <- "n"
    ticks <- c(.0001,.001,.01,.1, 1)
    p1 <- log(.025)
    p2 <- log(.01)
  }else{
    texty <- pstack + .03
    top <- 1
    logger <-  ""
    axer <- "s"
    p1 <- .025
    p2 <- .01
  }
  
#plot data
#this says which quadrant of the figure to use base on pos (1 to 4)
  co <- c(0,.5,0,.5)
  posit <- cbind(c(0,0,.5,.5), rep(.5,4),rep(0,4),c(.5,.5,0,0))
  co <- co+posit[,pos]
#Squishes the plots for the figure that also has the distrubution plots
  n <- 0
  if(short == TRUE){
    co <- co - c(0,0,0,0.075)
    n <- .015
    if(log == TRUE){
      texty <- texty*1.013  
    }else{
      texty <- texty+.013
    }
  }
  par(mar=c(3.5,4,2,1), fig = co, new = T, mgp = c(1.4,.4,0))
  
  plot(nodups, pstack, xlab = "ln(Repertoire Size)", ylab = "pValues",
       pch = 20, log = logger, yaxt = axer, ylim = c(.0001, top), cex.axis = .75,
       panel.first={
         abline(p1, 0, col = "grey50", lwd = 2, lty = 3)
         abline(p2, 0, col = "red", lwd = 2, lty = 3)
         })
  lines(nodups, pstack, col="black")
  title(titles, line = .15, cex.main = 1.1)
#adds the significance lines
  if(log == TRUE){
    axis(2, ticks, labels = ticks, cex.axis = .75)
  }
#plots the small group N
  text(nodups, texty+n, ntext, cex = .8)
  x <- which(co == 3)
  
#INSET
#this gives the inset that shows the significant p-vals on
#a log10 scale.  Otherwise, this code is more or less the same as above
  if(INSET == TRUE){
    sig <- which(pstack <= .05)
    if(length(sig)==0){
      stop("No significant p values.  Inset not plotted.")
    }
    highBou <- max(sig)
    lowBou <- min(sig)
    co2 <- c(.2, .37, .3, .45)
    co2 <- co2+posit[,pos]
    if(short == TRUE){
      co2 <- co2 - c(0,0,0.075,.075)  
    }
    par(mar=c(0,0,2.1,0), fig = co2, new = T)
    texty <- pstack*(2.5+n)
    top <- .25
    logger <- "y"
    axer <- "n"
    ticks <- c(.0001,.001,.01,.1)
    plot(nodups[lowBou:highBou], pstack[lowBou:highBou],
         xlab = "Repertoire Size", ylab = "pValues",
         pch = 20, log = logger, yaxt = axer, ylim = c(.0001, top),
         cex.axis = .6, cex= .8, xlim = c(nodups[lowBou], nodups[highBou]+.05),
         panel.first={
           abline(log10(.025), 0, col = "grey50", lwd = 2, lty = 3)
           abline(log10(.01), 0, col = "red", lwd = 2, lty = 3)
         })
    title("pValues < .05", line = .2, cex.main = .8)
    par(mgp = c(1.4,.6,0))
    axis(2, ticks, labels = ticks, cex.axis = .6, las = 2)
    text(nodups[lowBou:highBou], texty[lowBou:highBou], ntext[lowBou:highBou], cex = .6)
  }
}
PredInter <- function(dataset, Reper = NULL){
#pull the data subsets
  if(is.null(Reper) == FALSE){
    Set1 <- dataset[which(dataset$repsize < Reper),]
    Set2 <- dataset[which(dataset$repsize >= Reper),]
    pinames <- c("Population", "Small", "Large")
    OC <- FALSE
  }else{
    Set1 <- dataset[which(dataset$o.c == 1),]
    Set2 <- dataset[which(dataset$o.c == 2),]
    if(length(Set1$nbirds) < 3 || length(Set2$nbirds) < 3){
      stop("\n Not enough Birds for Prediction Intervals. \n 3 or more required! ")}
    pinames <- c("Population", "Song-Stable", "Song-Plastic")
    OC <- TRUE
  }
#run the three metanalyses
  PopulationPrin <- QuickMeta(dataset, Reper = Reper, ReturnALL = TRUE, OC = OC, pred = TRUE)
  OC <- FALSE
  Set1Prin <- QuickMeta(Set1, Reper = Reper, Popu = TRUE, ReturnALL = TRUE, OC = OC, pred = TRUE)
  Set2Prin <- QuickMeta(Set2, Reper = Reper, Popu = TRUE, ReturnALL = TRUE, OC = OC, pred = TRUE)
#pull confidence intervals
  P <- Datapuller(PopulationPrin, type = "PI")
  S1 <- Datapuller(Set1Prin, type = "PI")
  S2 <- Datapuller(Set2Prin, type = "PI")
  PIList<- list(P, S1, S2)
  names(PIList) <- pinames
  return(PIList)
}
Datapuller <- function(Meta, type = "p", Rep = FALSE){
  MetaPrint <- capture.output(print(Meta))
#What it says on the tin: this pulls data out of the
#meta-analysis printout that you can't otherwise access
#logic is the same for each section: 1) Find the line your
#data is on. 2) Pull the last elements(s), because that is
#where the numbers are. 3) clean the strings to make
#numeric if necessary, else just output the string
  
#pull p values for between groups
  if(type == "p"){
    line <- Linepuller ("Between groups", MetaPrint)
    lastelement <- length(line[[1]])
    Output <- as.numeric(line[[1]][lastelement])
#pull prediction intervals
  }else if(type == "PI"){
    pos <- grep("Prediction interval", MetaPrint)
    line <- strsplit(MetaPrint[pos], "\\s+")
    lastelement <- length(line[[1]])
    Output <- line[[1]][(lastelement-1):lastelement]
#remove garbage and make numeric
    Output <- gsub(paste0(";"), "",  Output)
    Output <- gsub(paste0("\\["), "",  Output)
    Output <- gsub(paste0("]"), "",  Output)
    Output <- as.numeric(Output)
#pull tau2 and I2
  }else if(type == "diag"){
    if(is.null(Rep) == TRUE){
      line <- Linepuller ("Group = Population", MetaPrint)
      lastelement <- length(line[[1]])
      I <- line[[1]][lastelement]
      Tau <- line[[1]][lastelement-1]
      Output <- list(c(I, Tau))
      return(Output)
    }else if(Rep == TRUE){
      line1 <- Linepuller ("Group = Smaller Syl Rep", MetaPrint)
      line2 <- Linepuller ("Group = Larger Syl Rep", MetaPrint)
      gnames <- c("Smaller", "Larger")
    }else{
      line1 <- Linepuller ("Group = Song-Stable", MetaPrint)
      line2 <- Linepuller ("Group = Song-Plastic", MetaPrint)
      gnames <- c("Song-Stable", "Song-Plastic")
    }
    lastelement1 <- length(line1[[1]])
    lastelement2 <- length(line2[[1]])
    I1 <- line1[[1]][lastelement1]
    Tau1 <- line1[[1]][lastelement1-1]
    I2 <- line2[[1]][lastelement2]
    Tau2 <- line2[[1]][lastelement2-1]
    Output <- list(c(I1, Tau1), c(I2, Tau2))
    names(Output) <- gnames
#Get the study weights
  }else if(type == "weight"){
    endrow <- length(Meta$studlab)+1
    lines <- strsplit(MetaPrint[2:endrow], "\\s+")
    index <- cumsum(sapply(lines, length))
    lines <- unlist(lines)
#the following -3 is based on the length of the group name (x Syl Rep = 3 chunks)
    Output <- as.numeric(lines[index-3])
#get the Z and associated p-val
  }else if(type == "signi"){
    line1 <- Linepuller ("Random effects model", MetaPrint)
    lastelement1 <- length(line1[[1]])
    if(length(line1) == 2){
      lastelement2 <- length(line1[[2]])
      pVal <- line1[[2]][lastelement2]
      Zval <- line1[[1]][lastelement1]
    }else{
      pVal <- line1[[1]][lastelement1]
      if(line1[[1]][lastelement1-1] == "<"){
        Zval <- line1[[1]][lastelement1-2]
      }else{Zval <- line1[[1]][lastelement1-1]}
    }
    Output <- c(Zval, pVal)
  }else{
    stop("\n Invalid type. \n type =") #'p', 'PI', 'diag', 'weight', or 'signi'.")
  }
  return(Output)
}
Linepuller <- function(string = "NULL", Meta){
#A subfunction for datapuller; splits the line into strings
  pos <- grep(paste(string), Meta)
  line <- strsplit(Meta[pos], "\\s+")
  return(line)
}
Ncounter <- function(fullrep, boundryrep){
#counts N for small rep group only!
  nreps <- length(fullrep)
#find the smallest and largest n
  smallestn <- length(which(fullrep < min(boundryrep)))
  largestn <-length(which(fullrep < max(boundryrep)))
#set up stuff for the loop
  ntext <- as.vector(NULL)
  counter <- smallestn-1
  adder <- 0
  adj <- smallestn - 1
#if1 deals with when you are at the last element
#if2 handles when a species is not duplicated
#else handles when a species is duplicated
  for(i in smallestn:largestn){
    counter <- counter+1
    if(counter == smallestn){
      ntext <- c(ntext, counter)
    }else if(boundryrep[i-adj] != boundryrep[i-smallestn]){
      counter <- counter + adder
      ntext <- c(ntext, counter)
      adder <- 0
    }else{
      counter <- counter-1
      adder <- adder+1
    }
  }
  return(ntext)
}
RepCutter <- function(rep, n){
#for the p-value plots:
#Makes sure you always have at least N studies in a group
#It is an issue, because some species have mutiple studies
#and those studies must pass the threshold as a unit
  rep <- sort(rep)
  last <- length(rep)
#Picks your smallest threshold
  if(rep[n] != rep[n+1]){
    rep <- rep[(n+1):last]
  }else{
    ind <- which(rep == rep[n])
    end <- max(ind)
    rep <- rep[(end+1):last]
  }
    last <- length(rep)
#picks your largest threshold
      if(rep[last-(n-1)] != rep[last-(n)]){
    cutrep <- rep[1:(last-(n-1))]
  }else{
    ind <- which(rep == rep[last-(n-1)])
    end <- min(ind)
    cutrep <- rep[1:(end)]
  }
  cutrep <- cutrep[which(duplicated(cutrep) == FALSE)]
  return(cutrep)
}
OCJackknife <- function(dataset, Spec = TRUE, Sort = TRUE){
#changes the OC grouping of a study or species
#sort data
  if(Sort == TRUE){
    dataset <- dataset[order(dataset$repsize),]
  }
#starter stuff for the loop 
  Species <- unique(dataset$species.name)
  pholder <- 1:length(Species)
#makes a vector where all of the OC groups are switched
  Switch <- sapply(dataset$o.c, Switcher)
#Species jackknife
  if(Spec == TRUE){
#create a new matrix where each col is the OC data
#with the grouping of one species swicthed
    Jack <- cbind(replicate(length(Species), dataset$o.c))
    JackKnife <- Replacer2(dataset, Switch, Jack)
    for (i in 1:length(Species)){
      dataset$o.c <- JackKnife[,i]
      pholder[i] <- QuickMeta(dataset, OC = TRUE, ReturnP = TRUE)
    }
    snames <- as.character(Species)
#Study jackknife
  }else if(Spec == FALSE){
#create a new matrix where each col is the OC data
#with the grouping of one study swicthed
    Jack <- cbind(replicate(length(dataset$o.c), dataset$o.c))
    JackKnife <- Replacer1(Jack, Switch)
    for (i in 1:length(dataset$o.c)){
      dataset$o.c <- JackKnife[,i]
      pholder[i] <- QuickMeta(dataset, OC = TRUE, ReturnP = TRUE)
    }
    snames <- as.character(dataset$citation)
  }
  Jacker <- cbind(snames, pholder)
  return(Jacker)
}
Switcher <- function(x){
#tests whether data is open(2) or closed(1) and creates a vector
#where everything is reversed
  if(x == 1){
    change <- 2
  }else{ change <- 1}
  return(change)
}
Replacer2 <- function(full, replace, final){
#Switches data tO->C or C->O in chunks based on species
  ind <- which(duplicated(full$species.name) == FALSE)
  counter <- 1
  adder <- 1
  j <- 1
#if1 handles when the last species in the category is a duplicate
#if2 handles when the current species is not a duplicate
#else handles when the current species is a duplicate
  for(i in 1:length(full$o.c)){
    if(i == length(full$o.c)){
      j <- j-1
      final[i,j] <-  replace[i]
    }else if(i == ind[j]){
      final[i,j] <-  replace[i]
      counter <- counter+adder
      adder <- 1
      j <- j+1
    }else{
      j <- j-1
      final[i,j] <-  replace[i]
      adder <- adder+1
      j <- j+1
    }
    
  }
  return(final)
}
Replacer1 <- function(full, replace){
#Switches data O->C or C->O study by study
  ind <- length(replace)
  for(i in 1:ind){
    full[i,i] <- replace[i]
  }
  return(full)
}


GroupSignificance <- function(dataset, Var = "Rep"){
#get the Z and p values for each subgroup and the full population 
#cut the data into subpopulations
  if(Var == "Rep"){
    group1 <- dataset[which(dataset$repsize < 38),]
    group2 <- dataset[which(dataset$repsize >= 38),]
    Reper <- 38
    OC <- FALSE
    gnames <- c("Population", "Small", "Large")
  }else if(Var == "OC"){
    group1 <- dataset[which(dataset$o.c == 1),]
    group2 <- dataset[which(dataset$o.c == 2),]
    OC <- TRUE
    Reper <- NULL
    gnames <- c("Population", "Song-Stable", "Song-Plastic")
  }
#run the meta-analyses
  p <- QuickMeta(dataset, Reper = Reper, OC = OC, ReturnALL = TRUE)
  OC <- FALSE
  g1 <- QuickMeta(group1, Reper = Reper, OC = OC, Popu = TRUE, ReturnALL = TRUE)
  g2 <- QuickMeta(group2, Reper = Reper, OC = OC, Popu = TRUE, ReturnALL = TRUE)
#pull the data form the printout
  Pop <- Datapuller(p, type = "signi")
  G1 <- Datapuller(g1, type = "signi")
  G2 <- Datapuller(g2, type = "signi")
#print it pretty
  Output <- list(Pop, G1, G2)
  names(Output) <- gnames
  return(Output)
}
PIStackedPlot <- function(dataset1, dataset2, dataset3){
#Pull prediction intervals
  oc1 <- unlist(PredInter(dataset1))
  LS1 <- unlist(PredInter(dataset1, Reper= 38))
  oc2 <- unlist(PredInter(dataset2))
  LS2 <- unlist(PredInter(dataset2, Reper= 38))
  oc3 <- unlist(PredInter(dataset3))
  LS3 <- unlist(PredInter(dataset3, Reper= 38))
  
#sort into line beginings and ends
  Starts <- c(oc1[c(1,3,5)], LS1[c(3,5)], oc2[c(1,3,5)],
              LS2[c(3,5)], oc3[c(1,3,5)], LS3[c(3,5)])
  Ends <- c(oc1[c(2,4,6)], LS1[c(4,6)], oc2[c(2,4,6)],
            LS2[c(4,6)], oc3[c(2,4,6)], LS3[c(4,6)])
#y position data and lables
  basepos <- c(1.1, 1.15, 1.2, 1.25, 1.3)
  allpos <- c(basepos, basepos-.5, basepos-1)
  Labels <- c("Population", "Stable", "Plastic", "Small", "Large")
  Labelpos <- c(rep(-.5, 5), rep(.5, 5), rep(0, 5))
  colors = c("black", "blue", "red", "cyan", "deeppink")
  
#Main plot
  par(mgp = c(1.5,.5,0))
  plot(0, type="n", ylab = "", xlab = "Correlation", xlim = c(-1, 1), yaxt = "n",
       ylim = (c(.05,1.45)), main="",
       panel.first = {
#Demarkation/grid lines
         segments(-1,1,1,1, lty = 'dashed', col = "grey68")
         segments(-1,.5,1,.5, lty = 'dashed', col = "grey68")
         segments(0,0,0,1.45, lty = 'dotted', col = "grey40")
         segments(.5,0,.5,1.45, lty = 'dotted', col = "grey68")
         segments(-.5,0,-.5,1.45, lty = 'dotted', col = "grey68")
         segments(-1,0,-1,1.45, lty = 'dotted', col = "grey85")
         segments(1,0,1,1.45, lty = 'dotted', col = "grey85")
       })
#Prediction interval lines
  segments(Starts, allpos, Ends, allpos, colors, lwd = 3)
#titles and edge markings ect
  text(Starts-.02,allpos,round(Starts, digits = 2), cex=.6, adj = 1)
  text(Ends+.02,allpos,round(Ends, digits = 2), cex=.6, adj = 0)
  text(0, c(1.4, .9, .4), c("Latency to Pairing Date", "Number of Offspring", "Number of Females"))
  text(-.9,1.4, "A", cex = 1, font = 2)
  text(-.9,.9, "B", cex = 1, font = 2)
  text(-.9,.4, "C", cex = 1, font = 2)
  text(-.35, 1.35, "<--Higher Mating Success", cex = .6, col = "grey40")
  text(.35, c(.85, .35), "Higher Mating Success-->", cex = .6, col = "grey40")
#Legend. :D
  legend(x = .75, y=1.4, legend=rev(Labels),
         col = rev(colors), pch = 15, xjust = 0, cex = .6)
}
StratifyTest <- function(dataset, reper = 38){
#get repertoire binarized
  dataset$repsize[which(dataset$repsize < reper)] <- 0
  dataset$repsize[which(dataset$repsize>= reper)] <- 1
  x1 <- dataset$repsize
  x2 <- dataset$o.c
  Output<- vector(mode="character", length=4) 
  
#stratify data and run t-tests
  Open <- dataset[which(x2 == 2),c(1,3,4,6)]
  Ol <- Open[which(Open$repsize == 1),4]
  Os <- Open[which(Open$repsize == 0),4]
  N2 <- c(length(Ol), length(Os))
  Otest <- t.test(Ol, Os)
  Output[1] <- Ttester(Otest, N2)
  
  Closed <- dataset[which(x2 == 1),c(1,3,4,6)]
  Cl <- Closed[which(Closed$repsize == 1),4]
  Cs <- Closed[which(Closed$repsize == 0),4]
  N2 <- c(length(Cl), length(Cs))
  Ctest <- t.test(Cl, Cs)
  Output[2] <- Ttester(Ctest, N2)
  
  
  Large <- dataset[which(x1 == 1),c(1,3,4,6)]
  Lc <- Large[which(Large$o.c == 1),4]
  Lo <- Large[which(Large$o.c == 2),4]
  N2 <- c(length(Lc), length(Lo))
  Ltest <- t.test(Lc, Lo)
  Output[3] <- Ttester(Ltest, N2)
  
  Small <- dataset[which(x1 == 0),c(1,3,4,6)]
  So <- Small[which(Small$o.c ==2),4]
  Sc <- Small[which(Small$o.c ==1),4]
  N2 <- c(length(Sc), length(So))
  Stest <- t.test(Sc, So)
  Output[4] <- Ttester(Stest, N2)
  return(Output)
}
Ttester <- function(n, n2){
#this just gives a pretty output for the t-tests in the stratify code
  paste(substr(n$data.name,1,2), "N:", n2[1], "mean:", round(n$estimate[1], digits = 4),
        substr(n$data.name,8,10), "N:", n2[2], "mean:", round(n$estimate[2], digits = 4),
        "T-stat:", round(n$statistic, digits = 4),
        "p-Value:", round(n$p.value, digits = 4), sep = " ")
}
ITAU <- function(dataset){
#get variable name
  dataname <- deparse(substitute(dataset))
#get the Tau2s and I2s
  datasetTauIRep38 <- QuickMeta(dataset, Reper = 38, Diag = TRUE)
  datasetIOC <- QuickMeta(dataset, OC = TRUE, Diag = TRUE)
  datasetTauIPOP <- QuickMeta(dataset, Popu = TRUE, Diag = TRUE)
#name the pop data
  names(datasetTauIPOP)<- paste(dataname,"Population", sep = " ")
#print it. :D
  print(c(datasetTauIPOP, datasetTauIRep38, datasetIOC))
}
WeightCompare <- function(dataset){
#get the what the average weight would be in fixed effect
  ave <- 100/length(dataset$citation)
#pull weights used and get range
  weightDS <- Datapuller(QuickMeta(dataset, Reper = 38, ReturnALL = TRUE), type="weight")
  rang <- range(weightDS)
  print(paste("Average:", ave, "  Range:", rang[1], "to", rang[2], sep = " "))
}
DistributionPlots <- function(PairDate, Females, Fledgling){
#pull data, combine, kill dups, sort, separate out closed-ended learners
  all <- rbind(PairDate, Females, Fledgling)
  nodups <-all[which(duplicated(all$species.name)==FALSE),]
  nodups <- nodups[order(nodups$repsize),]
  Cind <- which(nodups$o.c==1)
  Close <-nodups[Cind,]
#stuff for graphing
  ticks <- length(nodups$species.name)
  loc <- which(nodups$repsize == 38)-.5
  par(mfrow=c(1,2), mar = c(5,3.5,.5,1), mgp =c(2.2,.6,0), new = T, fig = c(0,.5,.5,1))
  legtop <- max(nodups$repsize) - .025*(max(nodups$repsize))

#plot for linear graph
  plot(0, type="n", xlim = c(1,ticks), ylim = c(min(nodups$repsize), max(nodups$repsize)),
       las = 2, xaxt = "n", xlab = "", ylab ="Repertoire Size", pch=21, col = "red",  cex = .6,
       cex.axis = .75,
       panel.first = {segments(loc, 0, loc, max(nodups$repsize), lty = 2, lwd = 1.5, col = "grey60")})
  segments(2:ticks-1, nodups$repsize[2:ticks-1], 2:ticks, nodups$repsize[2:ticks],
           lty = 1, lwd = 1.5, col = "black")
  points(1:ticks, nodups$repsize, pch=19, col = "red")
  points(Cind, Close$repsize, pch=19, col = "blue")
  axis(1, 1:ticks, labels = nodups$species.name, cex.axis = .6, las=2, font = 3)
  legend(1, legtop, legend=c("Song Stable", "Song Plastic"),
         col=c("blue", "red"), pch = 19, cex=.8)

#plot for ln() graph
  par(mar=c(0,0,2.1,0), mar = c(5,3.5,.5,1), mgp =c(2.2,.6,0), fig = c(.5, 1, .5, 1), new = T)
  nodups$repsize <- log(nodups$repsize)
  Close$repsize <- log(Close$repsize)
  plot(0, type="n", xlim = c(1,ticks), ylim = c(min(nodups$repsize), max(nodups$repsize)),
       las = 2, xaxt = "n", xlab = "", ylab ="ln(Repertoire Size)", pch=21, col = "red",  cex = .6,
       cex.axis = .75,
       panel.first = {segments(loc, 0, loc, max(nodups$repsize), lty = 2, lwd = 1.5, col = "grey60")})
  segments(2:ticks-1, nodups$repsize[2:ticks-1], 2:ticks, nodups$repsize[2:ticks],
           lty = 1, lwd = 1.5, col = "black")
  points(1:ticks, nodups$repsize, pch=19, col = "red")
  points(Cind, Close$repsize, pch=19, col = "blue")
  axis(1, 1:ticks, labels = nodups$species.name, cex.axis = .6, las=2, font  = 3)
}

