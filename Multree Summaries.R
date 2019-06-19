setwd("D:/Documents/R/2018-Meta-Analysis/Output/Rev3")
pop <- read.mulTree("Population")
songstab <- read.mulTree("SongStab")
sylcont <- read.mulTree("SylCont")
syldislow <- read.mulTree("SylDis 18.5")
syldismed <- read.mulTree("SylDis 38")
syldishigh <- read.mulTree("SylDis 216")

#hdr is the suggested method, but because SE is fixed at one, the 95% credint cannot
#be caculated for teh residual variance which makes the results for that weird.
#We went with teh quantile method since teh results were very similar, and the residual
#variance was appropriate estimated, even though we did not report it in ther paper.

#sink("Multree.txt")
#summary(pop, prob=95)
#summary(songstab, prob=95)
#summary(sylcont, prob=95)
#summary(syldislow, prob=95)
#summary(syldismed, prob=95)
#summary(syldishigh, prob=95)
#sink(NULL)

#sink("Multree.txt")
pop <- summary(pop, use.hdr=FALSE, prob=95)
songstab <- summary(songstab, use.hdr=FALSE, prob=95)
sylcont <- summary(sylcont, use.hdr=FALSE, prob=95)
syldislow <- summary(syldislow, use.hdr=FALSE, prob=95)
syldismed <- summary(syldismed, use.hdr=FALSE, prob=95)
syldishigh<- summary(syldishigh, use.hdr=FALSE, prob=95)
sink(NULL)

datable <- rbind(pop[1,], songstab[1:2,], sylcont[1:2,], syldislow[1:2,], syldismed[1:2,], syldishigh[1:2,])
rownames(datable) <- c("Population",
                          "Stable", "Plastic",
                          "Intercept", "Slope",
                          "Smaller < 18.5", "Larger >= 18.5",
                          "Smaller < 38", "Larger >= 38",
                          "Smaller < 216", "Larger >= 216")
datable <- round(datable, digits=3)
datable[,2:3]
datable[,2] <- apply(datable, 1, function(x) paste0("[", x[2], ";", x[3], "]"))
datable <- datable[,1:2]
colnames(datable) <- c("Post Mean", "95% CredInt")

xtable(datable)
