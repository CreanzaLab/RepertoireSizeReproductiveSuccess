########
#Coded by Cristina Robinson
#Last Modified 10-23-2017
#Written in R-Studio Version 1.1.383
#R Version 3.4.2
#phytools v0.6-30     ape v4.1    maps v3.2.0
########



rm(list = objects())
#install.packages("maps") 
#install.packages("ape")
#install.packages("phytools")
#Libraries
library(maps)
library(ape)
library(phytools)

#setwd()  ##set your directory here
dir <- getwd()
OLD <- FALSE #set this to TRUE if you are using a version of phytools older than v0.5-38
source(paste(dir, "/Tree functions_.R", sep = ""))

#Loading the trees takes a while, because the function is
#making a consensus tree from 1000 raw trees
setwd(paste(dir, "/Trees", sep =""))
Tree1 <- LoadPrettyTree("FledgeTree.tre")
Tree3 <- LoadPrettyTree("PairTree.tre")
#Tree2 <- LoadPrettyTree("FemaleTree.tre")
setwd(paste(dir))
Data1 <- LoadBirdData("Main Num Fledge Final.csv", Tree1)
Data3 <- LoadBirdData("Main Pair Date Final.csv", Tree3, Pair = TRUE)#reverses signs on pairdate correlations
#Data2 <- LoadBirdData("Main Num Females Final.csv", Tree2)

#set your desired output figure folder;
#the default is a (mostly) empty folder included in repository.
setwd(paste(dir, "/Output", sep = ""))

#Each plot takes a while, because it runs 5K sims per ANOVA and 2-3 ANOVAs
#the consol will print where you are at for those steps
set.seed(49)  ###This is my lucky number. :)
DoubleRainbowplot(Tree1, Data1[,1], Data1[,2], Data1[,3], Core = Data1[,4],
                  title="Number of Offspring", ANOVA = TRUE, Fig = "B")
DoubleRainbowplot(Tree3, Data3[,1], Data3[,2], Data3[,3], Core = Data3[,4],
                  title="Latency to Pairing Date", ANOVA = TRUE, Fig = "A")

#not included in final paper!!!
#DoubleRainbowplot(Tree2, Data2[,1], Data2[,2], Data2[,3], Core = Data2[,4],
#                  title="Number of Females")
