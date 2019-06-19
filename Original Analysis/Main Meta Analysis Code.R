########
#Coded by Cristina Robinson
#Last Modified 10-23-2017
#Written in R-Studio Version 1.1.383
#R Version 3.4.2
#metafor v2.0-0   Matrix v1.2-11   meta v4.8-4 
########

#As is, this code assumes you downloaded the full repository intact,
#and that you navigated to the folder containing the code.
#It pulls based on my folder organization scheme

rm(list = objects())
#install.packages("meta") 
#install.packages("matrix")
#install.packages("metafor")
library(meta)
library(Matrix)
library(metafor)

##set your directory here
#setwd()
dir <- getwd()
source(paste(dir, "/Main functions_.R", sep = ""))

#Read in primary data
PairDate <- read.csv("Main Pair Date Final.csv", header = TRUE)
Females <- read.csv("Main Num Females Final.csv", header = TRUE)
Offspring <- read.csv("Main Num Fledge Final.csv", header = TRUE)

#load territory controled dataset
#setwd(paste(dir, "/Territory Controlled Datasets", sep = ""))
#PairDate <- read.csv("Secondary Pair Date Final.csv", header = TRUE)
#Females <- read.csv("Secondary Num Females Final.csv", header = TRUE)
#Offspring <- read.csv("Secondary Num Fledge Final.csv", header = TRUE)
#dir <- getwd()

#set your desired output figure folder;
#the default is a (mostly) empty folder included in repository.
if (file.exists("Output")){
  setwd(file.path(dir, "Output"))
} else {
  dir.create(file.path(dir, "Output"))
  setwd(file.path(dir, "Output"))
  
}

#Forest Plots
QuickMeta(PairDate, Reper = 38, PrinPDF = TRUE)
QuickMeta(Females, Reper = 38, PrinPDF = TRUE)
QuickMeta(Offspring, Reper = 38, PrinPDF = TRUE)
QuickMeta(PairDate, OC = TRUE, PrinPDF = TRUE)
QuickMeta(Females, OC = TRUE, PrinPDF = TRUE)
QuickMeta(Offspring, OC = TRUE, PrinPDF = TRUE)

#Pulling the I^2 and Tau^2 (in that order!)
ITAU(PairDate)
ITAU(Offspring)
ITAU(Females)

#Compare study weights
WeightCompare(PairDate)
WeightCompare(Offspring)
WeightCompare(Females)

#Pulls the z value and p value (in that order!)
#for each population and subpopulation
#remember that the effectsizes and CI are plotted on the forest plots!
GroupSignificance(PairDate, Var = "Rep")
GroupSignificance(Females, Var = "Rep")
GroupSignificance(Offspring, Var = "Rep")
#These don't matter, because these groups were not significantly different
#GroupSignificance(PairDate, Var = "OC")
#GroupSignificance(Females, Var = "OC")
#GroupSignificance(Offspring, Var = "OC")


#Jackknife OC p values
###BY SPECIES
OCJackknife(PairDate)
OCJackknife(Offspring)
OCJackknife(Females)
###BY STUDY
OCJackknife(PairDate, Spec = FALSE, Sort = FALSE)
OCJackknife(Offspring, Spec = FALSE, Sort = FALSE)
OCJackknife(Females, Spec = FALSE, Sort = FALSE)


#Main pvalue and distrubution plots
Titles <-c("Latency to Pairing Date", "Number of Offspring")
pdf("Distribution and pValuePlot.pdf")
par(mar=c(0,0,0,0), fig=c(0,1,0,1), cex = .8, mgp = c(1.4,.4,0))
plot.new()
text(0,1, "A", cex = 1.5, font = 2)
text(.53,1, "B", cex = 1.5, font = 2)
text(0,.4, "C", cex = 1.5, font = 2)
text(.53,.4, "D", cex = 1.5, font = 2)
DistributionPlots(PairDate, Females, Offspring)
pvalueplot(PairDate, smallestn = 2, log = FALSE, INSET = TRUE, pos = 3, short = TRUE, title = Titles[1])
pvalueplot(Offspring, smallestn = 2, log = FALSE, INSET = TRUE, pos = 4, short = TRUE, title = Titles[2])
dev.off()

#Supplemental Pvalue Plots
#these are the datasets with the two most positive and negative correlations removed
LowestRemovedPD <- PairDate[order(PairDate$cor),][3:20,]
HighestRemovedPD <- PairDate[order(PairDate$cor),][1:18,]
LowestRemovedOF <- Offspring[order(Offspring$cor),][3:14,]
HighestRemovedOF <- Offspring[order(Offspring$cor),][1:12,]

Titles <-c("Latency to Pairing Date\nHighest Removed", "Latency to Pairing Date\nLowest Removed",
           "Number of Offspring\nHighest Removed", "Number of Offspring\nLowest Removed",
           "Number of Females") 
pdf("pValuePlot.pdf")
par(mar=c(0,0,0,0), fig=c(0,1,0,1), cex = .8, mgp = c(1.4,.4,0))
plot.new()
text(0,1, "A", cex = 1.5, font = 2)
text(.55,1, "B", cex = 1.5, font = 2)
text(0,.45, "C", cex = 1.5, font = 2)
text(.55,.45, "D", cex = 1.5, font = 2)
pvalueplot(HighestRemovedPD, smallestn = 2, log = FALSE, INSET = TRUE, pos = 1, title = Titles[1])
pvalueplot(LowestRemovedPD, smallestn = 2, log = FALSE, INSET = TRUE, pos = 2, title = Titles[2])
pvalueplot(HighestRemovedOF, smallestn = 2, log = FALSE, INSET = TRUE, pos = 3, title = Titles[3])
pvalueplot(LowestRemovedOF, smallestn = 2, log = FALSE, INSET = TRUE, pos = 4, title = Titles[4])
#Page 2
plot.new()
par(mar=c(0,0,0,0), fig=c(0,1,0,1))
pvalueplot(Females, smallestn = 2, log = FALSE, pos = 1, title = Titles[5])
dev.off()
#return to normalcy
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0))


#Stacked Prediction Interval Plot
pdf("StackedPredictionIntervals.pdf")
PIStackedPlot(PairDate, Offspring, Females)
dev.off()


#Stratify data to look for interactions
StratifyTest(PairDate)
StratifyTest(Offspring)
