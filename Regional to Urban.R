###################################################################################
#
# Ohio bee data
#
# Step 1: Ohio as Regional Pool to Urban Cleveland Pool
#
# CWM and Functional Diversity - Null Models
#
# KI Perry; 20 July 2021
# CA Shepard: 26 September 2022 
#CA Shepard: 20 February 2023 (Updated)
#CA Shepard: 3 October 2023 (update )
#CA Shepard: December 20th Update, based on code changes that occurred in a different version of the same document earlier in the year
###################################################################################

t <- read.csv("btraits_23.csv", row.names=1)
a <- read.csv("bcomm_25.localanalysis.csv", row.names=1)


str(a)
a1 <- a #save the original dataset
a <- a[2:362]
str(a)

rowSums(a1[2:362])
rowSums(a) #all sites have at least 10 species? 
#sites that have less than 10: t1-BE, T1C, t1H, t1SV, t6c, t7-c, t7-sv, t8-c 
#Sites that have 4 or fewer species: t1-BE,  t1H, t6c, t7-c
#The above 4 sites are removed below
a<-a[-60,]
a<-a[-57,]
a<-a[-54,]
a<-a[-49,]

rowSums(a)#all sites have 4+ species

#Need to check Meg. concinna and Meg. pusilla presence
#This is a taxonomy name change - all concinna in North America are actually pusilla
a$Megachileconcinna #Okay one concinna found
a$Megachilepusilla # No pusilla found

#This may cause a change in the datasets because M. pusilla is native, M. concinna is not

a$Megachilepusilla<- ifelse(a$Megachileconcinna != 0 | a$Megachilepusilla != 0, 1, 0)
a$Megachilepusilla #Okay yes, the above code made M. pusilla have a 1 in locations where either it or M concinna had a 1 

a<- a[,-239] #removing M. concinna


str(t)
names(t)
colnames(t) <- c("bl", "lec", "nest", "soc", "ori")
names(t)

plot(t)
cor(t, method = c("pearson"), use = "complete.obs")
str(t)

t1 <- t #save the original dataset

t$lec <- as.factor(t$lec)
t$nest <- as.factor(t$nest)
t$soc <- as.factor(t$soc)

str(t) # have to keep origin as a integer for the trait distance matrix to work

levels(t$lec)
levels(t$nest)
levels(t$soc)

#check body length for normality
hist(t$bl)
hist(log(t$bl))

t2 <- t #create another duplicate dataset before we transform
t2$bl <- log(t2$bl + 1)

t2<- t2[-239,] #removing M. concinna

#Double check that all species are present in both datasets
#Double check if a species is present in one dataset but not the other
setdiff(colnames(a), rownames(t2))
setdiff(rownames (t2), colnames(a))

rownames(t2) == colnames(a) # we are good to go!


#Going to go ahead and create CSV files for the regional pool with the correction to species name
#AND with the removed sites so that doesn't need to happen for each local site

write.csv(t2, file = "reg.traits.26.csv")

a2<-a1[-60,]
a2<-a2[-57,]
a2<-a2[-54,]
a2<-a2[-49,]
a$trmt<-a2$X.1
write.csv(a, file = "reg.comm.26.csv")

##############################################################################
## Observed Community Metrics----

if (!suppressWarnings(require(FD))) install.packages("FD")
citation("FD")

if (!suppressWarnings(require(picante))) install.packages("picante")
citation("picante")

if (!suppressWarnings(require(gawdis))) install.packages("gawdis")
citation("gawdis")

if (!suppressWarnings(require(betapart))) install.packages("betapart")
citation("betapart")


#observed CWM
cwm.obs <- functcomp(t2, as.matrix(a), CWM.type = "all")
cwm.obs

#observed taxonomic beta diversity
# create beta part object for analyses
str(a)
b.core <- betapart.core(a)

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
b.dist <- beta.pair(b.core, index.family = "sorensen")
str(b.dist)

#observed functional beta diversity
#create distance matrix with the traits
#optimized feature helps to weight the traits equally
tdis <- gawdis(t2, w.type = "optimized", opti.maxiter = 500)
attr(tdis, "correls")
attr(tdis, "weights")

# save trait weights for the null model
wt <- c(0.38, 0.15, 0.11, 0.10, 0.26)

#now run a principal coordinates analysis (PCoA) so we can collapse these traits into 
#a few continuous axes for the functional diversity calculations
pcoB <- dudi.hillsmith(as.matrix(tdis), scannf = FALSE, nf = 4)
pcoB

# check correlations among axes and traits
cor(pcoB$li, t1, use = "complete.obs")

sum(pcoB$eig[1:4]) / sum(pcoB$eig)
sum(pcoB$eig[1:3]) / sum(pcoB$eig)
sum(pcoB$eig[1:2]) / sum(pcoB$eig)

sum(pcoB$eig[2:3]) / sum(pcoB$eig)

t.ax <- as.matrix(pcoB$li[1:3])
b.fun <- functional.beta.pair(a, t.ax, index.family = "sorensen")
str(b.fun)

#We have now calculated all the indices with our observed bee data
#next, we need to run the null model!

#test the code to make sure it is working correctly
#single iteration of the null model
#Randomize community data matrix (maintains site species richness but not species frequency)
#Because only fixing species richness, susceptible to type 1 errors

#for (i in 1:nrow(a)){
 #randomizedB <- randomizeMatrix(samp = a, null.model = "richness", iterations = 1)
#}

#rowSums(a) == rowSums(randomizedB) #should be TRUE
#colSums(a) == colSums(randomizedB) #should be FALSE
#I think a few were true when there were 0 species
#colSums(a)
#colSums(randomizedB) #it worked!

#null CWMs
#cwm.test <- functcomp(t2, as.matrix(randomizedB), CWM.type = "all")
#cwm.test

#test randomizing the trait matrix as well
#traitsRand <- t2[sample(1:nrow(t2)),]
#rownames(traitsRand) <- rownames(t2)

#run the null models with 999 iterations
numberReps <- 999


#create matrices to store the results of each iteration of the null model, for each trait and index:
# for cwms
nbl <- nlec_0 <- nlec_1 <- nlec_2 <- nnest_1 <- nnest_2 <- nnest_3 <- nnest_4 <- nnest_5 <- nsoc_1 <- nsoc_2 <- nsoc_3 <- nsoc_4 <- nori_0 <- nori_1 <- matrix(NA,
      nrow = nrow(a), ncol = numberReps, dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# for taxonomic beta diversity
nbsim <- nbsne <- nbsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                  dimnames = list(rownames(a), paste0("n", 1:numberReps)))

# for functional beta diversity
nfsim <- nfsne <- nfsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                  dimnames = list(rownames(a), paste0("n", 1:numberReps)))

#create null model for each repetition:

for(i in 1:numberReps){
  print(i) 
  
  # randomized trait matrix
  ntraits <- t2[sample(1:nrow(t2)),]
  rownames(ntraits) <- rownames(t2)
  
  # randomized presence/absence matrix
  nsp <- randomizeMatrix(samp = a, null.model = "richness")
  
  # randomized trait distance matrix
  ntdis <- gawdis(ntraits, w.type = "user", W = wt)
  
  # CWM calculations
  cwm.null <- functcomp(x = ntraits, a = as.matrix(nsp), CWM.type = "all")
  nbl[,i] <- cwm.null$bl
  nlec_0[,i] <- cwm.null$lec_0
  nlec_1[,i] <- cwm.null$lec_1
  nlec_2[,i] <- cwm.null$lec_2
  nnest_1[,i] <- cwm.null$nest_1
  nnest_2[,i] <- cwm.null$nest_2
  nnest_3[,i] <- cwm.null$nest_3
  nnest_4[,i] <- cwm.null$nest_4
  nnest_5[,i] <- cwm.null$nest_5
  nsoc_1[,i] <- cwm.null$soc_1
  nsoc_2[,i] <- cwm.null$soc_2
  nsoc_3[,i] <- cwm.null$soc_3
  nsoc_4[,i] <- cwm.null$soc_4
  nori_0[,i] <- cwm.null$ori_0
  nori_1[,i] <- cwm.null$ori_1
  
  # Taxonomic beta diversity indices
  nb.core <- betapart.core(nsp)
  nb.dist <- beta.pair(nb.core, index.family = "sorensen")
  nsim.dist <- as.matrix(nb.dist$beta.sim)
  nsne.dist <- as.matrix(nb.dist$beta.sne)
  nsor.dist <- as.matrix(nb.dist$beta.sor)
  nbsim[,i] <- colMeans(nsim.dist)
  nbsne[,i] <- colMeans(nsne.dist)
  nbsor[,i] <- colMeans(nsor.dist)
  
  # Functional beta diversity indices
  npco <- dudi.hillsmith(as.matrix(ntdis), scannf = FALSE, nf = 3)
  nt <- as.matrix(npco$li)
  nb.fun <- functional.beta.pair(nsp, nt, index.family = "sorensen")
  nfsim.dist <- as.matrix(nb.fun$funct.beta.sim)
  nfsne.dist <- as.matrix(nb.fun$funct.beta.sne)
  nfsor.dist <- as.matrix(nb.fun$funct.beta.sor)
  nfsim[,i] <- colMeans(nfsim.dist)
  nfsne[,i] <- colMeans(nfsne.dist)
  nfsor[,i] <- colMeans(nfsor.dist)
  
}

# save the output matrices
write.csv(nbl, file = "Final.25.Regional to Urban_Nulls/nbl.csv")
write.csv(nlec_0, file = "Final.25.Regional to Urban_Nulls/nlec_0.csv")
write.csv(nlec_1, file = "Final.25.Regional to Urban_Nulls/nlec_1.csv")
write.csv(nlec_2, file = "Final.25.Regional to Urban_Nulls/nlec_2.csv")
write.csv(nnest_1, file = "Final.25.Regional to Urban_Nulls/nnest_1.csv")
write.csv(nnest_2, file = "Final.25.Regional to Urban_Nulls/nnest_2.csv")
write.csv(nnest_3, file = "Final.25.Regional to Urban_Nulls/nnest_3.csv")
write.csv(nnest_4, file = "Final.25.Regional to Urban_Nulls/nnest_4.csv")
write.csv(nnest_5, file = "Final.25.Regional to Urban_Nulls/nnest_5.csv")
write.csv(nsoc_1, file = "Final.25.Regional to Urban_Nulls/nsoc_1.csv")
write.csv(nsoc_2, file = "Final.25.Regional to Urban_Nulls/nsoc_2.csv")
write.csv(nsoc_3, file = "Final.25.Regional to Urban_Nulls/nsoc_3.csv")
write.csv(nsoc_4, file = "Final.25.Regional to Urban_Nulls/nsoc_4.csv")
write.csv(nori_0, file = "Final.25.Regional to Urban_Nulls/nori_0.csv")
write.csv(nori_1, file = "Final.25.Regional to Urban_Nulls/nori_1.csv")

write.csv(nbsim, file = "Final.25.Regional to Urban_Nulls/tbeta_sim.csv")
write.csv(nbsne, file = "Final.25.Regional to Urban_Nulls/tbeta_sne.csv")
write.csv(nbsor, file = "Final.25.Regional to Urban_Nulls/tbeta_sor.csv")

write.csv(nfsim, file = "Final.25.Regional to Urban_Nulls/fbeta_sim.csv")
write.csv(nfsne, file = "Final.25.Regional to Urban_Nulls/fbeta_sne.csv")
write.csv(nfsor, file = "Final.25.Regional to Urban_Nulls/fbeta_sor.csv")

# load the output matrices

nbl <- read.csv("Final.25.Regional to Urban_Nulls/nbl.csv", row.names=1)
nlec_0 <- read.csv("Final.25.Regional to Urban_Nulls/nlec_0.csv", row.names=1)
nlec_1 <- read.csv("Final.25.Regional to Urban_Nulls/nlec_1.csv", row.names=1)
nlec_2 <- read.csv("Final.25.Regional to Urban_Nulls/nlec_2.csv", row.names=1)
nnest_1 <- read.csv("Final.25.Regional to Urban_Nulls/nnest_1.csv", row.names=1)
nnest_2 <- read.csv("Final.25.Regional to Urban_Nulls/nnest_2.csv", row.names=1)
nnest_3 <- read.csv("Final.25.Regional to Urban_Nulls/nnest_3.csv", row.names=1)
nnest_4 <- read.csv("Final.25.Regional to Urban_Nulls/nnest_4.csv", row.names=1)
nnest_5 <- read.csv("Final.25.Regional to Urban_Nulls/nnest_5.csv", row.names=1)
nsoc_1 <- read.csv("Final.25.Regional to Urban_Nulls/nsoc_1.csv", row.names=1)
nsoc_2 <- read.csv("Final.25.Regional to Urban_Nulls/nsoc_2.csv", row.names=1)
nsoc_3 <- read.csv("Final.25.Regional to Urban_Nulls/nsoc_3.csv", row.names=1)
nsoc_4 <- read.csv("Final.25.Regional to Urban_Nulls/nsoc_4.csv", row.names=1)
nori_0 <- read.csv("Final.25.Regional to Urban_Nulls/nori_0.csv", row.names=1)
nori_1 <- read.csv("Final.25.Regional to Urban_Nulls/nori_1.csv", row.names=1)

nbsim <- read.csv("Final.25.Regional to Urban_Nulls/tbeta_sim.csv", row.names=1)
nbsne <- read.csv("Final.25.Regional to Urban_Nulls/tbeta_sne.csv", row.names=1)
nbsor <- read.csv("Final.25.Regional to Urban_Nulls/tbeta_sim.csv", row.names=1)

nfsim <- read.csv("Final.25.Regional to Urban_Nulls/fbeta_sim.csv", row.names=1)
nfsne <- read.csv("Final.25.Regional to Urban_Nulls/fbeta_sne.csv", row.names=1)
nfsor <- read.csv("Final.25.Regional to Urban_Nulls/fbeta_sor.csv", row.names=1)

## SES Calculations
#calculate standardized effect sizes (SES) for each trait and index
#the effect size is the difference between the observed value and the expected one
#then divide the effect size by the standard deviation of the null distribution to get the standardized effect size
#allows comparison among sites with different numbers of species

#test different methods with body length
#1
#meanNull_bl <- rowMeans(nbl)
#ES_bl <- cwm.obs$bl - meanNull_bl
#sdNull_bl <- apply(nbl, 1, sd)
#SES_bl.1 <- ES_bl / sdNull_bl

#2
#SES_bl.2 <-(cwm.obs$bl - rowMeans(nbl)) / apply(nbl, 1, sd, na.rm=T)

#3
#SES_bl.3 <-(cwm.obs$bl - apply(nbl, MARGIN = 1, mean)) / apply(nbl, MARGIN = 1, sd, na.rm=T)

#data.frame(SES_bl.1, SES_bl.2, SES_bl.3) #all three give the same results!

#################################################################################################################

# Community weighted means
cwm.obs.m <- as.matrix(t(colMeans(cwm.obs)))

 ## body length----
SES_bl <- (cwm.obs$bl - apply(nbl, MARGIN = 1, mean)) / apply(nbl, MARGIN = 1, sd, na.rm=T)
SES_bl

boxplot(SES_bl)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_bl, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nbl))
abline(v = cwm.obs.m[,1], col = "blue", lwd = 2)

pval.bl <- apply(cbind(cwm.obs$bl, nbl), MARGIN = 1, rank)[1,] / 1000
pval.bl
bl <- cbind(cwm.obs$bl, nbl)
bl.m <- as.matrix(t(colMeans(bl)))
pval.bl.a <- apply(bl.m, MARGIN = 1, rank)[1,] / 1000
pval.bl.a

# let's double check the p value calculationg above with a Wilcoxon Rank Test
w.bl <- wilcox.test(SES_bl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bl


## lec_0 - Kleptoparasitic----
SES_lec_0 <- (cwm.obs$lec_0 - apply(nlec_0, MARGIN = 1, mean)) / apply(nlec_0, MARGIN = 1, sd, na.rm=T)
SES_lec_0

boxplot(SES_lec_0, ylim = c(-4.0, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_lec_0, pch = 19, cex = 1.5, ylim = c(-4.0, 0.5))
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nlec_0))
abline (v = cwm.obs.m[,2], col = "blue", lwd = 2)

pval.lec_0 <- apply(cbind(cwm.obs$lec_0, nlec_0), MARGIN = 1, rank)[1,] / 1000
pval.lec_0
lec_0 <- cbind(cwm.obs$lec_0, nlec_0)
lec_0.m <- as.matrix(t(colMeans(lec_0)))
pval.lec_0.a <- apply(lec_0.m, MARGIN = 1, rank)[1,] / 1000
pval.lec_0.a

w.lec_0 <- wilcox.test(SES_lec_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.lec_0


## lec_1 - Generalist----
SES_lec_1 <- (cwm.obs$lec_1 - apply(nlec_1, MARGIN = 1, mean)) / apply(nlec_1, MARGIN = 1, sd, na.rm=T)
SES_lec_1

boxplot(SES_lec_1, ylim = c(-0.5, 6.0))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_lec_1, pch = 19, cex = 1.5, ylim = c(-0.5, 6.0))
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nlec_1))
abline (v = cwm.obs.m[,3], col = "blue", lwd = 2)

pval.lec_1 <- apply(cbind(cwm.obs$lec_1, nlec_1), MARGIN = 1, rank)[1,] / 1000
pval.lec_1
lec_1 <- cbind(cwm.obs$lec_1, nlec_1)
lec_1.m <- as.matrix(t(colMeans(lec_1)))
pval.lec_1.a <- apply(lec_1.m, MARGIN = 1, rank)[1,] / 1000
pval.lec_1.a

w.lec_1 <- wilcox.test(SES_lec_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.lec_1


## lec_2 - Specialist----
SES_lec_2 <- (cwm.obs$lec_2 - apply(nlec_2, MARGIN = 1, mean)) / apply(nlec_2, MARGIN = 1, sd, na.rm=T)
SES_lec_2

boxplot(SES_lec_2, ylim = c(-3.5, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_lec_2, pch = 19, cex = 1.5, ylim = c(-3.5, 0.5))
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nlec_2))
abline (v = cwm.obs.m[,4], col = "blue", lwd = 2)

pval.lec_2 <- apply(cbind(cwm.obs$lec_2, nlec_2), MARGIN = 1, rank)[1,] / 1000
pval.lec_2
lec_2 <- cbind(cwm.obs$lec_2, nlec_2)
lec_2.m <- as.matrix(t(colMeans(lec_2)))
pval.lec_2.a <- apply(lec_2.m, MARGIN = 1, rank)[1,] / 1000
pval.lec_2.a

w.lec_2 <- wilcox.test(SES_lec_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.lec_2


## nest_1 - Soil----
SES_nest_1 <- (cwm.obs$nest_1 - apply(nnest_1, MARGIN = 1, mean)) / apply(nnest_1, MARGIN = 1, sd, na.rm=T)
SES_nest_1

boxplot(SES_nest_1)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_nest_1, pch = 19, cex = 1.5)
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nnest_1))
abline (v = cwm.obs.m[,5], col = "blue", lwd = 2)

pval.nest_1 <- apply(cbind(cwm.obs$nest_1, nnest_1), MARGIN = 1, rank)[1,] / 1000
pval.nest_1
nest_1 <- cbind(cwm.obs$nest_1, nnest_1)
nest_1.m <- as.matrix(t(colMeans(nest_1)))
pval.nest_1.a <- apply(nest_1.m, MARGIN = 1, rank)[1,] / 1000
pval.nest_1.a

w.nest_1 <- wilcox.test(SES_nest_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest_1


## nest_2 - Cavity----
SES_nest_2 <- (cwm.obs$nest_2 - apply(nnest_2, MARGIN = 1, mean)) / apply(nnest_2, MARGIN = 1, sd, na.rm=T)
SES_nest_2

boxplot(SES_nest_2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_nest_2, pch = 19, cex = 1.5)
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nnest_2))
abline (v = cwm.obs.m[,6], col = "blue", lwd = 2)

pval.nest_2 <- apply(cbind(cwm.obs$nest_2, nnest_2), MARGIN = 1, rank)[1,] / 1000
pval.nest_2
nest_2 <- cbind(cwm.obs$nest_2, nnest_2)
nest_2.m <- as.matrix(t(colMeans(nest_2)))
pval.nest_2.a <- apply(nest_2.m, MARGIN = 1, rank)[1,] / 1000
pval.nest_2.a

w.nest_2 <- wilcox.test(SES_nest_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest_2


## nest_3 - Colony----
SES_nest_3 <- (cwm.obs$nest_3 - apply(nnest_3, MARGIN = 1, mean)) / apply(nnest_3, MARGIN = 1, sd, na.rm=T)
SES_nest_3

boxplot(SES_nest_3)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_nest_3, pch = 19, cex = 1.5)
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nnest_3))
abline (v = cwm.obs.m[,7], col = "blue", lwd = 2)

pval.nest_3 <- apply(cbind(cwm.obs$nest_3, nnest_3), MARGIN = 1, rank)[1,] / 1000
pval.nest_3
nest_3 <- cbind(cwm.obs$nest_3, nnest_3)
nest_3.m <- as.matrix(t(colMeans(nest_3)))
pval.nest_3.a <- apply(nest_3.m, MARGIN = 1, rank)[1,] / 1000
pval.nest_3.a

w.nest_3 <- wilcox.test(SES_nest_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest_3


## nest_4 - Pithy Stems----
SES_nest_4 <- (cwm.obs$nest_4 - apply(nnest_4, MARGIN = 1, mean)) / apply(nnest_4, MARGIN = 1, sd, na.rm=T)
SES_nest_4

boxplot(SES_nest_4)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_nest_4, pch = 19, cex = 1.5)
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nnest_4))
abline (v = cwm.obs.m[,8], col = "blue", lwd = 2)

pval.nest_4 <- apply(cbind(cwm.obs$nest_4, nnest_4), MARGIN = 1, rank)[1,] / 1000
pval.nest_4
nest_4 <- cbind(cwm.obs$nest_4, nnest_4)
nest_4.m <- as.matrix(t(colMeans(nest_4)))
pval.nest_4.a <- apply(nest_4.m, MARGIN = 1, rank)[1,] / 1000
pval.nest_4.a

w.nest_4 <- wilcox.test(SES_nest_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest_4



## nest_5 - Wood----
SES_nest_5 <- (cwm.obs$nest_5 - apply(nnest_5, MARGIN = 1, mean)) / apply(nnest_5, MARGIN = 1, sd, na.rm=T)
SES_nest_5

boxplot(SES_nest_5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_nest_5, pch = 19, cex = 1.5)
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nnest_5))
abline (v = cwm.obs.m[,9], col = "blue", lwd = 2)

pval.nest_5 <- apply(cbind(cwm.obs$nest_5, nnest_5), MARGIN = 1, rank)[1,] / 1000
pval.nest_5
nest_5 <- cbind(cwm.obs$nest_5, nnest_5)
nest_5.m <- as.matrix(t(colMeans(nest_5)))
pval.nest_5.a <- apply(nest_5.m, MARGIN = 1, rank)[1,] / 1000
pval.nest_5.a

w.nest_5 <- wilcox.test(SES_nest_5, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.nest_5


## soc_1 - Subsocial----
SES_soc_1 <- (cwm.obs$soc_1 - apply(nsoc_1, MARGIN = 1, mean)) / apply(nsoc_1, MARGIN = 1, sd, na.rm=T)
SES_soc_1

boxplot(SES_soc_1)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_soc_1, pch = 19, cex = 1.5)
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nsoc_1))
abline (v = cwm.obs.m[,10], col = "blue", lwd = 2)

pval.soc_1 <- apply(cbind(cwm.obs$soc_1, nsoc_1), MARGIN = 1, rank)[1,] / 1000
pval.soc_1
soc_1 <- cbind(cwm.obs$soc_1, nsoc_1)
soc_1.m <- as.matrix(t(colMeans(soc_1)))
pval.soc_1.a <- apply(soc_1.m, MARGIN = 1, rank)[1,] / 1000
pval.soc_1.a

w.soc_1 <- wilcox.test(SES_soc_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.soc_1


## soc_2 - Solitary----
SES_soc_2 <- (cwm.obs$soc_2 - apply(nsoc_2, MARGIN = 1, mean)) / apply(nsoc_2, MARGIN = 1, sd, na.rm=T)
SES_soc_2

boxplot(SES_soc_2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_soc_2, pch = 19, cex = 1.5)
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nsoc_2))
abline (v = cwm.obs.m[,11], col = "blue", lwd = 2)

pval.soc_2 <- apply(cbind(cwm.obs$soc_2, nsoc_2), MARGIN = 1, rank)[1,] / 1000
pval.soc_2
soc_2 <- cbind(cwm.obs$soc_2, nsoc_2)
soc_2.m <- as.matrix(t(colMeans(soc_2)))
pval.soc_2.a <- apply(soc_2.m, MARGIN = 1, rank)[1,] / 1000
pval.soc_2.a

w.soc_2 <- wilcox.test(SES_soc_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.soc_2


## soc_3 - Eusocial----
SES_soc_3 <- (cwm.obs$soc_3 - apply(nsoc_3, MARGIN = 1, mean)) / apply(nsoc_3, MARGIN = 1, sd, na.rm=T)
SES_soc_3

boxplot(SES_soc_3)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_soc_3, pch = 19, cex = 1.5)
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nsoc_3))
abline (v = cwm.obs.m[,12], col = "blue", lwd = 2)

pval.soc_3 <- apply(cbind(cwm.obs$soc_3, nsoc_3), MARGIN = 1, rank)[1,] / 1000
pval.soc_3
soc_3 <- cbind(cwm.obs$soc_3, nsoc_3)
soc_3.m <- as.matrix(t(colMeans(soc_3)))
pval.soc_3.a <- apply(soc_3.m, MARGIN = 1, rank)[1,] / 1000
pval.soc_3.a

w.soc_3 <- wilcox.test(SES_soc_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.soc_3


## soc_4 - Parasitic----
SES_soc_4 <- (cwm.obs$soc_4 - apply(nsoc_4, MARGIN = 1, mean)) / apply(nsoc_4, MARGIN = 1, sd, na.rm=T)
SES_soc_4

boxplot(SES_soc_4, ylim = c(-4.0, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_soc_4, pch = 19, cex = 1.5, ylim = c(-4.0, 0.5))
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nsoc_4))
abline (v = cwm.obs.m[,13], col = "blue", lwd = 2)

pval.soc_4 <- apply(cbind(cwm.obs$soc_4, nsoc_4), MARGIN = 1, rank)[1,] / 1000
pval.soc_4
soc_4 <- cbind(cwm.obs$soc_4, nsoc_4)
soc_4.m <- as.matrix(t(colMeans(soc_4)))
pval.soc_4.a <- apply(soc_4.m, MARGIN = 1, rank)[1,] / 1000
pval.soc_4.a

w.soc_4 <- wilcox.test(SES_soc_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.soc_4


## ori_0 - Native----
SES_ori_0 <- (cwm.obs$ori_0 - apply(nori_0, MARGIN = 1, mean)) / apply(nori_0, MARGIN = 1, sd, na.rm=T)
SES_ori_0

boxplot(SES_ori_0, ylim = c(-9.0, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_ori_0, pch = 19, cex = 1.5, ylim = c(-9.0, 0.5))
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nori_0))
abline (v = cwm.obs.m[,14], col = "blue", lwd = 2)

pval.ori_0 <- apply(cbind(cwm.obs$ori_0, nori_0), MARGIN = 1, rank)[1,] / 1000
pval.ori_0
ori_0 <- cbind(cwm.obs$ori_0, nori_0)
ori_0.m <- as.matrix(t(colMeans(ori_0)))
pval.ori_0.a <- apply(ori_0.m, MARGIN = 1, rank)[1,] / 1000
pval.ori_0.a

w.ori_0 <- wilcox.test(SES_ori_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.ori_0


## ori_2 - Exotic----
SES_ori_1 <- (cwm.obs$ori_1 - apply(nori_1, MARGIN = 1, mean)) / apply(nori_1, MARGIN = 1, sd, na.rm=T)
SES_ori_1

boxplot(SES_ori_1, ylim = c(-0.5, 9.0))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_ori_1, pch = 19, cex = 1.5, ylim = c(-0.5, 9.0))
abline (h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nori_1))
abline (v = cwm.obs.m[,15], col = "blue", lwd = 2)

pval.ori_1 <- apply(cbind(cwm.obs$ori_1, nori_1), MARGIN = 1, rank)[1,] / 1000
pval.ori_1
ori_1 <- cbind(cwm.obs$ori_1, nori_1)
ori_1.m <- as.matrix(t(colMeans(ori_1)))
pval.ori_1.a <- apply(ori_1.m, MARGIN = 1, rank)[1,] / 1000
pval.ori_1.a

w.ori_1 <- wilcox.test(SES_ori_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.ori_1


###############################################
## Diversity Indices##

# Taxonomic diversity----

beta.sor <- as.matrix(b.dist$beta.sor)
beta.sor <- colMeans(beta.sor)

beta.sim <- as.matrix(b.dist$beta.sim)
beta.sim <- colMeans(beta.sim)

beta.sne <- as.matrix(b.dist$beta.sne)
beta.sne <- colMeans(beta.sne)

beta.t <- data.frame(beta.sor, beta.sim, beta.sne)
beta.t.m <- as.matrix(t(colMeans(beta.t)))

## taxonomic diveristy - beta sor----
SES_bsor <- (beta.t$beta.sor - apply(nbsor, MARGIN = 1, mean)) / apply(nbsor, MARGIN = 1, sd, na.rm=T)
SES_bsor

boxplot(SES_bsor)
boxplot(SES_bsor, ylim=c(-75, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_bsor, pch = 19, cex = 1.5, ylim=c(-75, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nbsor), main = "Distribution of expected values - Taxonomic Beta-Diversity", xlim=c(0.5, 1))
abline(v = beta.t.m[,1], col = "blue", lwd = 2)

pval.beta.sor <- apply(cbind(beta.t$beta.sor, nbsor), MARGIN = 1, rank)[1,] / 1000
pval.beta.sor
tbsor <- cbind(beta.t$beta.sor, nbsor)
tbsor.m <- as.matrix(t(colMeans(tbsor)))
pval.beta.sor.a <- apply(tbsor.m, MARGIN = 1, rank)[1,] / 1000
pval.beta.sor.a

w.bsor <- wilcox.test(SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsor


## taxonomic diversity - beta sim----
SES_bsim <- (beta.t$beta.sim - apply(nbsim, MARGIN = 1, mean)) / apply(nbsim, MARGIN = 1, sd, na.rm=T)
SES_bsim

boxplot(SES_bsim)
boxplot(SES_bsim, ylim = c(-75, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_bsim, pch = 19, cex = 1.5, ylim = c(-75, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nbsim), main = "Distribution of expected values - Taxonomic Beta-Diversity - Turnover", xlim=c(0.4, 1))
abline(v = beta.t.m[,2], col = "blue", lwd = 2)

pval.beta.sim <- apply(cbind(beta.t$beta.sim, nbsim), MARGIN = 1, rank)[1,] / 1000
pval.beta.sim
tbsim <- cbind(beta.t$beta.sim, nbsim)
tbsim.m <- as.matrix(t(colMeans(tbsim)))
pval.beta.sim.a <- apply(tbsim.m, MARGIN = 1, rank)[1,] / 1000
pval.beta.sim.a

w.bsim <- wilcox.test(SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsim


## taxonomic diversity - beta sne----
SES_bsne <- (beta.t$beta.sne - apply(nbsne, MARGIN = 1, mean)) / apply(nbsne, MARGIN = 1, sd, na.rm=T)
SES_bsne

boxplot(SES_bsne)
boxplot(SES_bsne, ylim=c(-0.5, 60))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_bsne, pch = 19, cex = 1.5, ylim=c(-0.5, 60))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nbsne), main = "Distribution of expected values - Taxonomic Beta-Diversity - Nestedness", xlim=c(0, 0.2))
abline(v = beta.t.m[,3], col = "blue", lwd = 2)

pval.beta.sne <- apply(cbind(beta.t$beta.sne, nbsne), MARGIN = 1, rank)[1,] / 1000
pval.beta.sne
tbsne <- cbind(beta.t$beta.sne, nbsne)
tbsne.m <- as.matrix(t(colMeans(tbsne)))
pval.beta.sne.a <- apply(tbsne.m, MARGIN = 1, rank)[1,] / 1000
pval.beta.sne.a

w.bsne <- wilcox.test(SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.bsne


# Functional Diversity----
fbeta.sor <- as.matrix(b.fun$funct.beta.sor)
fbeta.sor <- colMeans(fbeta.sor)

fbeta.sim <- as.matrix(b.fun$funct.beta.sim)
fbeta.sim <- colMeans(fbeta.sim)

fbeta.sne <- as.matrix(b.fun$funct.beta.sne)
fbeta.sne <- colMeans(fbeta.sne)

beta.f <- data.frame(fbeta.sor, fbeta.sim, fbeta.sne)
beta.f.m <- as.matrix(t(colMeans(beta.f)))

plot(fbeta.sor, pch = 19, cex = 1.5)
plot(fbeta.sim, pch = 19, cex = 1.5)
plot(fbeta.sne, pch = 19, cex = 1.5)

## functional diversity - beta sor-----
SES_fbsor <- (beta.f$fbeta.sor - apply(nfsor, MARGIN = 1, mean)) / apply(nfsor, MARGIN = 1, sd, na.rm=T)
SES_fbsor

boxplot(SES_fbsor, ylim = c(-0.5, 10))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_fbsor, pch = 19, cex = 1.5, ylim = c(-0.5, 10))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nfsor), main = "Distribution of expected values - Functional Beta-Diversity", xlim = c(0, 1))
abline(v = beta.f.m[,1], col = "blue", lwd = 2)

pval.fbeta.sor <- apply(cbind(beta.f$fbeta.sor, nfsor), MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sor
fbsor <- cbind(beta.f$fbeta.sor, nfsor)
fbsor.m <- as.matrix(t(colMeans(fbsor)))
pval.fbeta.sor.a <- apply(fbsor.m, MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sor.a

w.fbsor <- wilcox.test(SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsor


## functional diversity - beta sim----
SES_fbsim <- (beta.f$fbeta.sim - apply(nfsim, MARGIN = 1, mean)) / apply(nfsim, MARGIN = 1, sd, na.rm=T)
SES_fbsim

boxplot(SES_fbsim)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_fbsim, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nfsim), main = "Distribution of expected values - Functional Beta-Diversity - Turnover", xlim = c(0, 0.3))
abline(v = beta.f.m[,2], col = "blue", lwd = 2)

pval.fbeta.sim <- apply(cbind(beta.f$fbeta.sim, nfsim), MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sim
fbsim <- cbind(beta.f$fbeta.sim, nfsim)
fbsim.m <- as.matrix(t(colMeans(fbsim)))
pval.fbeta.sim.a <- apply(fbsim.m, MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sim.a

w.fbsim <- wilcox.test(SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsim


## functional diversity - beta sne----
SES_fbsne <- (beta.f$fbeta.sne - apply(nfsne, MARGIN = 1, mean)) / apply(nfsne, MARGIN = 1, sd, na.rm=T)
SES_fbsne

boxplot(SES_fbsne)
boxplot(SES_fbsne, ylim = c(-0.5, 8))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

plot(SES_fbsne, pch = 19, cex = 1.5, ylim = c(-0.5, 8))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

hist(as.matrix(nfsne), main = "Distribution of expected values - Functional Beta-Diversity - Turnover")
abline(v = beta.f.m[,3], col = "blue", lwd = 2)

pval.fbeta.sne <- apply(cbind(beta.f$fbeta.sne, nfsne), MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sne
fbsne <- cbind(beta.f$fbeta.sne, nfsne)
fbsne.m <- as.matrix(t(colMeans(fbsne)))
pval.fbeta.sne.a <- apply(fbsne.m, MARGIN = 1, rank)[1,] / 1000
pval.fbeta.sne.a

w.fbsne <- wilcox.test(SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
w.fbsne


###################################################################################
# combine the SES values from each trait type into a matrix so we can graph the results----

if (!suppressWarnings(require(viridis))) install.packages("viridis")
citation("viridis")

if (!suppressWarnings(require(reshape2))) install.packages("reshape2")
citation("reshape2")

SES_traits <- as.data.frame(cbind(SES_ori_1, SES_ori_0, SES_lec_2, SES_lec_1, SES_lec_0, SES_bl))
colnames(SES_traits) <- c("Exotic", "Native", "Specialist", "Generalist", "Kleptoparasitic" , "Body Length")
SES_traits <- melt(SES_traits)
colnames(SES_traits) <- c("trait","ses")

par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = SES_traits, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Traits")
stripchart(ses ~ trait, data = SES_traits, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)



## Nesting Traits

SES_Nest <- as.data.frame(cbind(SES_nest_5, SES_nest_4, SES_nest_3, SES_nest_2, SES_nest_1))
colnames(SES_Nest) <- c("Wood","Pithy Stems", "Hive", "Cavity", "Soil")
SES_Nest <- melt(SES_Nest)
colnames(SES_Nest) <- c("nest","ses")

par(mar = c(5,9,4,2))
boxplot(ses ~ nest, data = SES_Nest, col = viridis(5, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Nesting Traits by Resource")
stripchart(ses ~ nest, data = SES_Nest, col = viridis(5),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)



## Sociality Traits
SES_Soc <- as.data.frame(cbind(SES_soc_4, SES_soc_3, SES_soc_1, SES_soc_2))
colnames(SES_Soc) <- c("Parasitic","Eusocial", "Subsocial", "Solitary")
SES_Soc <- melt(SES_Soc)
colnames(SES_Soc) <- c("soc","ses")

par(mar = c(5,9,4,2))
boxplot(ses ~ soc, data = SES_Soc, col = viridis(4, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Sociality Traits")
stripchart(ses ~ soc, data = SES_Soc, col = viridis(4),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)


## Taxonomic Beta-Diversity
SES_TBeta <- as.data.frame(cbind(SES_bsne, SES_bsim, SES_bsor))
colnames(SES_TBeta) <- c("Nestedness", "Turnover", "Total Beta-Diversity")
SES_TBeta <- melt(SES_TBeta)
colnames(SES_TBeta) <- c("tbeta","ses")

par(mar = c(5,10,4,2))
boxplot(ses ~ tbeta, data = SES_TBeta, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Taxonomic Beta-Diversity")
stripchart(ses ~ tbeta, data = SES_TBeta, col = viridis(3),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)


## Functional Beta-Diversity
SES_FBeta <- as.data.frame(cbind(SES_fbsne, SES_fbsim, SES_fbsor))
colnames(SES_FBeta) <- c("Nestedness","Turnover", "Total Beta-Diversity")
SES_FBeta <- melt(SES_FBeta)
colnames(SES_FBeta) <- c("fbeta","ses")

par(mar = c(5,10,4,2))
boxplot(ses ~ fbeta, data = SES_FBeta, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Functional Beta-Diversity")
stripchart(ses ~ fbeta, data = SES_FBeta, col = viridis(3),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)



## figure panel - Diversity Indices

### figure panel - all diversity indices

library(ggplot2)
#install.packages("ggthemes")
library(ggthemes)

SES_TBeta$metric <- rep(c("Taxonomic"),each = 189)
SES_FBeta$metric <- rep(c("Functional"),each = 189)

names(SES_TBeta) <- names(SES_FBeta) 
#The above code was added to get rid of the "names do not match" error
SES_div <- rbind(SES_TBeta, SES_FBeta)
SES_div.m <- melt(SES_div)
colnames(SES_div.m) <- c("metric", "var", "ses")

SES_div.m$metric <- factor(SES_div.m$metric, levels = c("Taxonomic", "Functional"))


#Don't use this figure, use the next one down
png("SES_Div.png", width = 2000, height = 1000, pointsize = 20)

ggplot(SES_div.m, aes(x=ses, y=var, fill = var)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(position = position_jitter(width = 0.2, height = 0.2),
               dotsize = 1,
               binaxis = "y",
               stackdir = "center") +
  facet_grid(metric ~ .) + 
  coord_cartesian(xlim = c(-100, 100)) +
  theme_few() +
  theme(text = element_text(size = 24, color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(vjust = 1, size = 28),
        axis.title.y = element_text(vjust = 1, size = 28),
        strip.text = element_text(face = "bold", size = rel(1.1)),
        strip.background = element_rect(fill = "gray88", color = "black"),
        legend.position = "none") +
  labs(x = "Standardized Effect Sizes (SES)", y = "Beta-diversity Metrics") +
  geom_vline(xintercept = 0, size = 1.2, linetype = "dashed") +
  scale_fill_viridis(alpha = 0.7, discrete = TRUE, option = "D")

dev.off()


png("Figures/Fig4.Regionaldivpanel.png", width = 3000, height = 1000, pointsize = 20)

par(mfrow=c(1,2)) # indicates one row, two columns
par(mar=c(5,15,4,2))

boxplot(ses ~ fbeta, data = SES_TBeta, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-90,80), cex.lab = 2, cex.axis = 1.45) #cex.main= 1.8(That is what I would use tomake the graph title size, but no titles))
stripchart(ses ~ fbeta, data = SES_TBeta, col = viridis(3),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#I am going to put in the text labels in adobe indesign since I can't get them placed consistently across all of the graphs


boxplot(ses ~ fbeta, data = SES_FBeta, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-3,8), cex.lab = 2, cex.axis = 1.45)
stripchart(ses ~ fbeta, data = SES_FBeta, col = viridis(3),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(6.5, 3.4, "B", pos = 4, font = 2, cex = 2.6)

dev.off()

#text(1700, 17.5, "A", pos = 4, font = 2, cex = 3)


## figure panel - Traits

png("Figures/SES_Traits.png", width = 600, height = 1600, pointsize = 20)

par(mfrow=c(3,1)) # indicates two rows, two columns
par(mar=c(5,10,4,2))

#traits
boxplot(ses ~ trait, data = SES_traits, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, cex.lab = 2, cex.axis = 1.7)
stripchart(ses ~ trait, data = SES_traits, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(6.2, 6, "A", pos = 4, font = 2, cex = 2.6)

#nesting
boxplot(ses ~ nest, data = SES_Nest, col = viridis(5, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, cex.lab = 2, cex.axis = 1.7)
stripchart(ses ~ nest, data = SES_Nest, col = viridis(5),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(4, 5, "B", pos = 4, font = 2, cex = 2.6)

#sociality
boxplot(ses ~ soc, data = SES_Soc, col = viridis(4, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, cex.lab = 2, cex.axis = 1.7)
stripchart(ses ~ soc, data = SES_Soc, col = viridis(4),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(4.5, 4.2, "C", pos = 4, font = 2, cex = 2.6)

dev.off()

#length and origin
SES_lengthandorigin <- as.data.frame(cbind( SES_ori_1, SES_ori_0, SES_bl ))
colnames(SES_lengthandorigin) <- c( "Exotic", "Native", "Body Length")
SES_lengthandorigin <- melt(SES_lengthandorigin)
colnames(SES_lengthandorigin) <- c("trait","ses")

par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = SES_lengthandorigin, col = viridis(5, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-8,8))
stripchart(ses ~ trait, data = SES_lengthandorigin, col = viridis(5),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

#lecty
SES_lect <- as.data.frame(cbind( SES_lec_2, SES_lec_1, SES_lec_0))
colnames(SES_lect) <- c( "Specialist", "Generalist", "Kleptoparasitic" )
SES_lect <- melt(SES_lect)
colnames(SES_lect) <- c("trait","ses")

par(mar = c(5,9,4,2))
boxplot(ses ~ trait, data = SES_lect, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-8,8))
stripchart(ses ~ trait, data = SES_lect, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)



## Nesting Traits

SES_Nest <- as.data.frame(cbind(SES_nest_5, SES_nest_4, SES_nest_3, SES_nest_2, SES_nest_1))
colnames(SES_Nest) <- c("Wood","Pithy Stems", "Hive", "Cavity", "Soil")
SES_Nest <- melt(SES_Nest)
colnames(SES_Nest) <- c("nest","ses")

par(mar = c(5,9,4,2))
boxplot(ses ~ nest, data = SES_Nest, col = viridis(5, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-8,8))
stripchart(ses ~ nest, data = SES_Nest, col = viridis(5),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)



## Sociality Traits
SES_Soc <- as.data.frame(cbind(SES_soc_4, SES_soc_3, SES_soc_1, SES_soc_2))
colnames(SES_Soc) <- c("Parasitic","Eusocial", "Subsocial", "Solitary")
SES_Soc <- melt(SES_Soc)
colnames(SES_Soc) <- c("soc","ses")

par(mar = c(5,9,4,2))
boxplot(ses ~ soc, data = SES_Soc, col = viridis(4, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-8,8))
stripchart(ses ~ soc, data = SES_Soc, col = viridis(4),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)



#Same graphs but Put in figure order for making paper figures----
png("Figures/Figure 5 Regionalsp traits panel.png", width = 1500, height = 1000, pointsize = 20)

par(mfrow=c(2,2)) # indicates two rows, two columns
par(mar = c(5,9,4,2)) # sets the margins around the figure

#length and origin
SES_lengthandorigin <- as.data.frame(cbind( SES_ori_1, SES_ori_0, SES_bl ))
colnames(SES_lengthandorigin) <- c( "Alien", "Native", "Body Length")
SES_lengthandorigin <- melt(SES_lengthandorigin)
colnames(SES_lengthandorigin) <- c("trait","ses")



boxplot(ses ~ trait, data = SES_lengthandorigin, col = viridis(5, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-9,9), cex.lab = 2, cex.axis = 1.45)
stripchart(ses ~ trait, data = SES_lengthandorigin, col = viridis(5),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(6.5, 3.3, "A", pos = 4, font = 2, cex = 2)


## Nesting Traits

SES_Nest <- as.data.frame(cbind(SES_nest_5, SES_nest_4, SES_nest_3, SES_nest_2, SES_nest_1))
colnames(SES_Nest) <- c("Wood","Pithy Stems", "Colony", "Cavity", "Soil")
SES_Nest <- melt(SES_Nest)
colnames(SES_Nest) <- c("nest","ses")

boxplot(ses ~ nest, data = SES_Nest, col = viridis(5, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-9,9), cex.lab = 2, cex.axis = 1.45)
stripchart(ses ~ nest, data = SES_Nest, col = viridis(5),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(6.5, 5.2, "B", pos = 4, font = 2, cex = 2)

#lecty
SES_lect <- as.data.frame(cbind( SES_lec_2, SES_lec_1, SES_lec_0))
colnames(SES_lect) <- c( "Specialist", "Generalist", "Kleptoparasitic" )
SES_lect <- melt(SES_lect)
colnames(SES_lect) <- c("trait","ses")


boxplot(ses ~ trait, data = SES_lect, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-9,9), cex.lab = 2, cex.axis = 1.45)
stripchart(ses ~ trait, data = SES_lect, col = viridis(6),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

#text(6.5, 3.3, "C", pos = 4, font = 2, cex = 2)


## Sociality Traits
SES_Soc <- as.data.frame(cbind(SES_soc_4, SES_soc_3, SES_soc_1, SES_soc_2))
colnames(SES_Soc) <- c("Parasitic","Eusocial", "Subsocial", "Solitary")
SES_Soc <- melt(SES_Soc)
colnames(SES_Soc) <- c("soc","ses")


boxplot(ses ~ soc, data = SES_Soc, col = viridis(4, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-9,9), cex.lab = 2, cex.axis = 1.45)
stripchart(ses ~ soc, data = SES_Soc, col = viridis(4),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(6.5, 4.2, "D", pos = 4, font = 2, cex = 2)

dev.off()



#Intro stats----
ordered(colSums(a))
sort(colSums(a), decreasing=T)
absent<-which(colSums(a)==0)
present<-a[,-absent]
#okay got that to work- 136 present species
sort(colSums(present), decreasing = T)


t$ori

Nonnative <- t[which(t$ori == 1),]
Nonnative



#Discussion stats
SES.all <- as.data.frame(cbind(SES_bl, SES_lec_0, SES_lec_1, SES_lec_2, SES_ori_0, SES_ori_1, 
                               SES_nest_1, SES_nest_2, SES_nest_3, SES_nest_4, SES_nest_5,
                               SES_soc_1, SES_soc_2, SES_soc_3, SES_soc_4, SES_bsor, SES_bsim,
                               SES_bsne, SES_fbsor, SES_fbsim, SES_fbsne))
write.csv(SES.all, file = "Final.25.Regional to Urban_Nulls/AllSESvals.csv")

#The below code seems left over from something
#nestsites <- read.csv("Highestfunctnestsites.csv", row.names=1)



#str(nestsites)
#nestsites1 <- nestsites #save the original dataset
#nestsites <- nestsites[2:362]
s#tr(nestsites)
#
#colSums(nestsites1[2:362])
#nestsites <- nestsites[,which(colSums(nestsites)!=0)] 
#colSums(nestsites) 
#Removed all species columsn where no member of the species was collected. 
#nestednessspp<-colnames(nestsites)
#nestednessspp

#nestednesstraits<-t2
#nestednesstraits = nestednesstraits[rownames(nestednesstraits) %in% nestednessspp, ]
#Got rid of all the species that were not part of the nestsites
#setdiff(colnames(nestsites), rownames(nestednesstraits))
#setdiff(rownames (nestednesstraits), colnames(nestsites))

#rownames(nestednesstraits) == colnames(nestsites) # we are good to go!

## pull out data for each treatment
#functnest <- nestsites[5:7,]
#functnest <- functnest[,which(colSums(functnest)!=0)] 
#colSums(functnest)
#functnestspp<-colnames(functnest)
#functnestednesstraits = nestednesstraits[rownames(nestednesstraits) %in% functnestspp, ]

#taxnest <- nestsites[2:4,]
#taxnest <- taxnest[,which(colSums(taxnest)!=0)] 
#colSums(taxnest)
#taxnestspp<-colnames(taxnest)
#taxnestednesstraits = nestednesstraits[rownames(nestednesstraits) %in% taxnestspp, ]
#summary(taxnestednesstraits)
#summary(functnestednesstraits)



##Figure for presentation----
par(mar = c(5,9,4,2)) # sets the margins around the figure
windows()
#length and origin
SES_lengthandorigin.1 <- as.data.frame(cbind( SES_ori_1, SES_ori_0 ))
colnames(SES_lengthandorigin.1) <- c( "Alien", "Native")
SES_lengthandorigin.1 <- melt(SES_lengthandorigin.1)
colnames(SES_lengthandorigin.1) <- c("trait","ses")



boxplot(ses ~ trait, data = SES_lengthandorigin.1, col = viridis(2, alpha = 0.55, begin = 0.4),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-9,9), cex.lab = 2, cex.axis = 1.45)
stripchart(ses ~ trait, data = SES_lengthandorigin.1, col = viridis(2, begin = 0.4),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(6.5, 3.3, "A", pos = 4, font = 2, cex = 2)


## Nesting Traits

SES_Nest.1 <- as.data.frame(cbind( SES_nest_2, SES_nest_1))
colnames(SES_Nest.1) <- c( "Cavity", "Soil")
SES_Nest.1 <- melt(SES_Nest.1)
colnames(SES_Nest.1) <- c("nest","ses")

boxplot(ses ~ nest, data = SES_Nest.1, col = viridis(3, alpha = 0.55, begin = 0),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-9,9), cex.lab = 2, cex.axis = 1.45)
stripchart(ses ~ nest, data = SES_Nest.1, col = viridis(3, begin = 0),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(6.5, 5.2, "B", pos = 4, font = 2, cex = 2)


## Taxonomic Beta-Diversity
#SES_TBeta <- as.data.frame(cbind(SES_bsne, SES_bsim, SES_bsor))
#colnames(SES_TBeta) <- c("Nestedness", "Turnover", "Total Beta-Diversity")
#SES_TBeta <- melt(SES_TBeta)
#colnames(SES_TBeta) <- c("tbeta","ses")
windows()
par(mar = c(5,10,4,2))
boxplot(SES_bsor, col = viridis(1, alpha = 0.55, begin = 0.75),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-60,60),main = "Taxonomic Beta-Diversity")
stripchart(SES_bsor, data = SES_TBeta, col = viridis(1, begin = 0.75),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)


## Functional Beta-Diversity
#SES_FBeta <- as.data.frame(cbind(SES_fbsne, SES_fbsim, SES_fbsor))
#colnames(SES_FBeta) <- c("Nestedness","Turnover", "Total Beta-Diversity")
#SES_FBeta <- melt(SES_FBeta)
#colnames(SES_FBeta) <- c("fbeta","ses")

par(mar = c(5,10,4,2))
boxplot(SES_fbsor, col = viridis(1, alpha = 0.55, begin = 0.8),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim= c(-6,7),main = "Functional Beta-Diversity")
stripchart(SES_fbsor, col = viridis(1, begin = 0.8),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
