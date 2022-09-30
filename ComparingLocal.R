###################################################################################
#
# Ohio bee data
#
# Step 2: Urban Cleveland Pool to Local Greenspaces
# Comparing among Local habitats
# 
#CWM and Functional Diversity - Null Models
# GLMMs
#
# KI Perry; 21 July 2021
#CA Shepard 07 June 2022
#CA Shepard; 26 September 2022
###################################################################################
#Creating the datasets----

t <- read.csv("./btraits_final.csv", row.names=1)
a <- read.csv("./bcomm_final.localanalysis.csv", row.names=1)


str(a)
a1 <- a #save the original dataset

# pull out treatments that we want to keep
farm <- a[which(a$trmt == "Farm"),]
str(farm)
Control <- a[which(a$trmt == "Control"),]
str(Control)
T1 <- a[which(a$trmt == "T1"),]
str(T1)
T8 <- a[which(a$trmt == "T8"),]
str(T8)

# create new dataset with only farm, Control, t1, and t8 treatments
a <- rbind(farm, Control,T1, T8)
str(a)

a <- a[2:360]
str(a)
rowSums(a1[2:360])

# create a vector with the column sums for each species
# species not found in Cleveland will have a 0
sp <- colSums(a)
sp

# removes any columns (i.e. species) that are not found in Cleveland
a <- a[, colSums(a != 0) > 0]
colSums(a)

str(t)
names(t)
colnames(t) <- c("bl", "lec", "nest", "soc", "ori")
names(t)

# add the sp vector as a column in the trait matrix, shows which species
# are found in Cleveland and which are absent (i.e. with a 0)
t$sp <- sp
t

# uses the sp values to remove rows of species not collected in Cleveland
# then remove the column because we don't need it anymore
t <- t[t$sp != 0, ]
t <- t[,-6]

plot(t)
cor(t, method = c("pearson"), use = "complete.obs")
str(t)

t1 <- t #save the original dataset

t$lec <- as.factor(t$lec)
t$nest <- as.factor(t$nest)
t$soc <- as.factor(t$soc)

str(t) # have to keep origin as a integer for the trait distance matrix to work

#check body length for normality
hist(t$bl)
hist(log(t$bl))

t2 <- t #create another duplicate dataset before we transform
t2$bl <- log(t2$bl + 1)
str(t2)

#Double check that all species are present in both datasets
#Double check if a species is present in one dataset but not the other
setdiff(colnames(a), rownames(t2))
setdiff(rownames (t2), colnames(a))

rownames(t2) == colnames(a) # we are good to go!

##############################################################################
## Observed Community Metrics----

#Loading needed packages ------
if (!suppressWarnings(require(FD))) install.packages("FD")
citation("FD")

if (!suppressWarnings(require(picante))) install.packages("picante")
citation("picante")

if (!suppressWarnings(require(gawdis))) install.packages("gawdis")
citation("gawdis")

if (!suppressWarnings(require(betapart))) install.packages("betapart")
citation("betapart")

rowSums(a)

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
wt <- c(0.38, 0.18, 0.11, 0.12, 0.20)

#now run a principal coordinates analysis (PCoA) so we can collapse these traits into 
#a few continuous axes for the functional diversity calculations
pcoB <- dudi.hillsmith(as.matrix(tdis), scannf = FALSE, nf = 4)
pcoB

# check correlations among axes and traits
cor(pcoB$li, t1, use = "complete.obs")

sum(pcoB$eig[1:4]) / sum(pcoB$eig)
sum(pcoB$eig[1:3]) / sum(pcoB$eig)
sum(pcoB$eig[1:2]) / sum(pcoB$eig)

t.ax <- as.matrix(pcoB$li[1:3])
b.fun <- functional.beta.pair(a, t.ax, index.family = "sorensen")
str(b.fun)

#observed functional alpha diversity
#run the rao function first!
bb.rao <- Rao(sample = t(a), dfunc = tdis, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
falpha <- bb.rao$FD$Alpha
falpha

#We have now calculated all the indices with our observed bee data
#next, we need to run the null model!----

#using independent swap method for randomizing the presence/absence matrix
#this will constrain the null communities by species richness and species frequency
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
nfalpha <- nfsim <- nfsne <- nfsor <- matrix(NA, nrow = nrow(a), ncol = numberReps, 
                                  dimnames = list(rownames(a), paste0("n", 1:numberReps)))

#create null model for each repetition:

for(i in 1:numberReps){
  print(i) 
  
  # randomized trait matrix
  ntraits <- t2[sample(1:nrow(t2)),]
  rownames(ntraits) <- rownames(t2)
  
  # randomized presence/absence matrix
  nsp <- randomizeMatrix(samp = a, null.model = "independentswap")
  
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
  
  # Functional alpha diversity
  nrao <- Rao(sample = t(nsp), dfunc = ntdis, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
  nfalpha[,i] <- nrao$FD$Alpha
  
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
write.csv(nbl, file = "nbl.csv")
write.csv(nlec_0, file = "nlec_0.csv")
write.csv(nlec_1, file = "nlec_1.csv")
write.csv(nlec_2, file = "nlec_2.csv")
write.csv(nnest_1, file = "nnest_1.csv")
write.csv(nnest_2, file = "nnest_2.csv")
write.csv(nnest_3, file = "nnest_3.csv")
write.csv(nnest_4, file = "nnest_4.csv")
write.csv(nnest_5, file = "nnest_5.csv")
write.csv(nsoc_1, file = "nsoc_1.csv")
write.csv(nsoc_2, file = "nsoc_2.csv")
write.csv(nsoc_3, file = "nsoc_3.csv")
write.csv(nsoc_4, file = "nsoc_4.csv")
write.csv(nori_0, file = "nori_0.csv")
write.csv(nori_1, file = "nori_1.csv")

write.csv(nbsim, file = "tbeta_sim.csv")
write.csv(nbsne, file = "tbeta_sne.csv")
write.csv(nbsor, file = "tbeta_sor.csv")

write.csv(nfalpha, file = "falpha.csv")

write.csv(nfsim, file = "fbeta_sim.csv")
write.csv(nfsne, file = "fbeta_sne.csv")
write.csv(nfsor, file = "fbeta_sor.csv")

# load the output matrices

nbl <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nbl.csv", row.names=1)
nlec_0 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nlec_0.csv", row.names=1)
nlec_1 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nlec_1.csv", row.names=1)
nlec_2 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nlec_2.csv", row.names=1)
nnest_1 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nnest_1.csv", row.names=1)
nnest_2 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nnest_2.csv", row.names=1)
nnest_3 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nnest_3.csv", row.names=1)
nnest_4 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nnest_4.csv", row.names=1)
nnest_5 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nnest_5.csv", row.names=1)
nsoc_1 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nsoc_1.csv", row.names=1)
nsoc_2 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nsoc_2.csv", row.names=1)
nsoc_3 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nsoc_3.csv", row.names=1)
nsoc_4 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nsoc_4.csv", row.names=1)
nori_0 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nori_0.csv", row.names=1)
nori_1 <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/nori_1.csv", row.names=1)

nfalpha <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/falpha.csv", row.names=1)

nbsim <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/tbeta_sim.csv", row.names=1)
nbsne <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/tbeta_sne.csv", row.names=1)
nbsor <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/tbeta_sim.csv", row.names=1)

nfsim <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/fbeta_sim.csv", row.names=1)
nfsne <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/fbeta_sne.csv", row.names=1)
nfsor <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes/fbeta_sor.csv", row.names=1)

# SES Calculations----
#calculate standardized effect sizes (SES) for each trait and index
#the effect size is the difference between the observed value and the expected one
#then divide the effect size by the standard deviation of the null distribution to get the standardized effect size
#allows comparison among sites with different numbers of species


## community weighted means----
## body length
SES_bl <- (cwm.obs$bl - apply(nbl, MARGIN = 1, mean)) / apply(nbl, MARGIN = 1, sd, na.rm=T)
SES_bl

## lec_0 - Kleptoparasitic
SES_lec_0 <- (cwm.obs$lec_0 - apply(nlec_0, MARGIN = 1, mean)) / apply(nlec_0, MARGIN = 1, sd, na.rm=T)
SES_lec_0

## lec_1 - Generalist
SES_lec_1 <- (cwm.obs$lec_1 - apply(nlec_1, MARGIN = 1, mean)) / apply(nlec_1, MARGIN = 1, sd, na.rm=T)
SES_lec_1

## lec_2 - Specialist
SES_lec_2 <- (cwm.obs$lec_2 - apply(nlec_2, MARGIN = 1, mean)) / apply(nlec_2, MARGIN = 1, sd, na.rm=T)
SES_lec_2

## nest_1 - Soil
SES_nest_1 <- (cwm.obs$nest_1 - apply(nnest_1, MARGIN = 1, mean)) / apply(nnest_1, MARGIN = 1, sd, na.rm=T)
SES_nest_1

## nest_2 - Cavity
SES_nest_2 <- (cwm.obs$nest_2 - apply(nnest_2, MARGIN = 1, mean)) / apply(nnest_2, MARGIN = 1, sd, na.rm=T)
SES_nest_2

## nest_3 - Hive
SES_nest_3 <- (cwm.obs$nest_3 - apply(nnest_3, MARGIN = 1, mean)) / apply(nnest_3, MARGIN = 1, sd, na.rm=T)
SES_nest_3

## nest_4 - Pithy Stems
SES_nest_4 <- (cwm.obs$nest_4 - apply(nnest_4, MARGIN = 1, mean)) / apply(nnest_4, MARGIN = 1, sd, na.rm=T)
SES_nest_4

## nest_5 - Wood
SES_nest_5 <- (cwm.obs$nest_5 - apply(nnest_5, MARGIN = 1, mean)) / apply(nnest_5, MARGIN = 1, sd, na.rm=T)
SES_nest_5

## soc_1 - Subsocial
SES_soc_1 <- (cwm.obs$soc_1 - apply(nsoc_1, MARGIN = 1, mean)) / apply(nsoc_1, MARGIN = 1, sd, na.rm=T)
SES_soc_1

## soc_2 - Solitary
SES_soc_2 <- (cwm.obs$soc_2 - apply(nsoc_2, MARGIN = 1, mean)) / apply(nsoc_2, MARGIN = 1, sd, na.rm=T)
SES_soc_2

## soc_3 - Eusocial
SES_soc_3 <- (cwm.obs$soc_3 - apply(nsoc_3, MARGIN = 1, mean)) / apply(nsoc_3, MARGIN = 1, sd, na.rm=T)
SES_soc_3

## soc_4 - Parasitic
SES_soc_4 <- (cwm.obs$soc_4 - apply(nsoc_4, MARGIN = 1, mean)) / apply(nsoc_4, MARGIN = 1, sd, na.rm=T)
SES_soc_4

## ori_0 - Native
SES_ori_0 <- (cwm.obs$ori_0 - apply(nori_0, MARGIN = 1, mean)) / apply(nori_0, MARGIN = 1, sd, na.rm=T)
SES_ori_0

## ori_2 - Exotic
SES_ori_1 <- (cwm.obs$ori_1 - apply(nori_1, MARGIN = 1, mean)) / apply(nori_1, MARGIN = 1, sd, na.rm=T)
SES_ori_1

## taxonomic beta diversity----
beta.sor <- as.matrix(b.dist$beta.sor)
beta.sor <- colMeans(beta.sor)

beta.sim <- as.matrix(b.dist$beta.sim)
beta.sim <- colMeans(beta.sim)

beta.sne <- as.matrix(b.dist$beta.sne)
beta.sne <- colMeans(beta.sne)

beta.t <- data.frame(beta.sor, beta.sim, beta.sne)

## taxonomic diveristy - beta sor
SES_bsor <- (beta.t$beta.sor - apply(nbsor, MARGIN = 1, mean)) / apply(nbsor, MARGIN = 1, sd, na.rm=T)
SES_bsor

## taxonomic diveristy - beta sim
SES_bsim <- (beta.t$beta.sim - apply(nbsim, MARGIN = 1, mean)) / apply(nbsim, MARGIN = 1, sd, na.rm=T)
SES_bsim

## taxonomic diveristy - beta sne
SES_bsne <- (beta.t$beta.sne - apply(nbsne, MARGIN = 1, mean)) / apply(nbsne, MARGIN = 1, sd, na.rm=T)
SES_bsne

## functional diversity----
fbeta.sor <- as.matrix(b.fun$funct.beta.sor)
fbeta.sor <- colMeans(fbeta.sor)

fbeta.sim <- as.matrix(b.fun$funct.beta.sim)
fbeta.sim <- colMeans(fbeta.sim)

fbeta.sne <- as.matrix(b.fun$funct.beta.sne)
fbeta.sne <- colMeans(fbeta.sne)

beta.f <- data.frame(falpha, fbeta.sor, fbeta.sim, fbeta.sne)

## functional alpha  - Rao
SES_falpha <- (beta.f$falpha - apply(nfalpha, MARGIN = 1, mean)) / apply(nfalpha, MARGIN = 1, sd, na.rm=T)
SES_falpha

## functional diveristy - beta sor
SES_fbsor <- (beta.f$fbeta.sor - apply(nfsor, MARGIN = 1, mean)) / apply(nfsor, MARGIN = 1, sd, na.rm=T)
SES_fbsor

## functional diveristy - beta sim
SES_fbsim <- (beta.f$fbeta.sim - apply(nfsim, MARGIN = 1, mean)) / apply(nfsim, MARGIN = 1, sd, na.rm=T)
SES_fbsim

## functional diveristy - beta sne
SES_fbsne <- (beta.f$fbeta.sne - apply(nfsne, MARGIN = 1, mean)) / apply(nfsne, MARGIN = 1, sd, na.rm=T)
SES_fbsne

## combine all indices into one matrix
SES.all <- as.data.frame(cbind(SES_bl, SES_lec_0, SES_lec_1, SES_lec_2, SES_ori_0, SES_ori_1, 
                           SES_nest_1, SES_nest_2, SES_nest_3, SES_nest_4, SES_nest_5,
                           SES_soc_1, SES_soc_2, SES_soc_3, SES_soc_4, SES_bsor, SES_bsim,
                           SES_bsne, SES_falpha, SES_fbsor, SES_fbsim, SES_fbsne))

write.csv(SES.all, file = "SES_Local.csv")
#import the SES data
SES.all <- read.csv("./SES_Local.csv", row.names = 1)

# import the landscape data-----
land <- read.csv("./landscape.localanalysis.csv", row.names = 1)
str(land)

land$trmt <- as.factor(land$trmt)
str(land)

## merge landscape data with SES data
SES <- merge(SES.all, land, by = c("row.names"))
str(SES)

SES$site <- as.factor(SES$Row.names)
str(SES)

## pull out data for each treatment
farm <- SES[which(SES$trmt == "Farm"),]
str(farm)

control <- SES[which(SES$trmt == "Control"),]
str(control)

T1 <- SES[which(SES$trmt == "T1"),]
str(T1)

T8 <- SES[which(SES$trmt == "T8"),]
str(T8)

#Loading more needed packages----
if (!suppressWarnings(require(nortest))) install.packages("nortest")
citation("nortest")

if (!suppressWarnings(require(car))) install.packages("car")
citation("car")

if (!suppressWarnings(require(emmeans))) install.packages("emmeans")
citation("emmeans")

if (!suppressWarnings(require(bbmle))) install.packages("bbmle")
citation("bbmle")

if (!suppressWarnings(require(DHARMa))) install.packages("DHARMa")
citation("DHARMa")

if (!suppressWarnings(require(lme4))) install.packages("lme4")
citation("lme4")

if (!suppressWarnings(require(ggplot2))) install.packages("ggplot2")
citation("ggplot2")

if (!suppressWarnings(require(sjPlot))) install.packages("sjPlot")
citation("sjPlot")

if (!suppressWarnings(require(jtools))) install.packages("jtools")
citation("jtools")

if (!suppressWarnings(require(interactions))) install.packages("interactions")
citation("interactions")


## check relationships among landscape variables----
dotchart(land$pland, pch = 19)
dotchart(land$lpi, pch = 19)
dotchart(land$enn, pch = 19)
plot(land[2:4], pch = 19)
cor(land[2:4], method = c("pearson"), use = "complete.obs")

# now let's compare by treatment----

## body length
hist(SES$SES_bl)
plot(SES$SES_bl, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bl ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
bl.farm <- wilcox.test(farm$SES_bl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bl.farm
bl.control <- wilcox.test(control$SES_bl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bl.control
bl.t1 <- wilcox.test(T1$SES_bl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bl.t1
bl.t8 <- wilcox.test(T8$SES_bl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bl.t8

## compare among treatments and landscape variables
dotchart(SES$SES_bl, group = SES$trmt, pch = 19)

##Need to write code for not comparing all at once----
FS_bl <- as.data.frame(cbind(farm$SES_bl, control$SES_bl))
#okay that makes a matrix, but it doesn't work because the different treatments have different locations. can't be compared equally
FS_comps<-SES[1:29,]
FS_comps<-FS_comps[-3,]
FS_comps<-FS_comps[-3,]
FS_comps<-FS_comps[-13,]
FS_comps<-FS_comps[-13,]
FS_comps<-FS_comps[-13,]
FS_comps<-FS_comps[-13,]

FS_comps<-FS_comps[-14,]
FS_comps<-FS_comps[-14,]
FS_comps<-FS_comps[-15,]
FS_comps<-FS_comps[-15,]
FS_comps<-FS_comps[-15,]
FS_comps<-FS_comps[-15,]
#Okay so the above created a matrix of Francis' results
with(FS_comps, bartlett.test(SES_bl ~ trmt))
#####IT WORKS!!!!!
#Okay so to do
#1- make a matrix of Katie results
#2- Make a matrix of both vacant lot results
#rewrite code for comparisons to be about each thing

#back to normal comparisons 
with(SES, bartlett.test(SES_bl ~ trmt))
with(SES, ad.test(SES_bl))

bl.mod.full <- glm(SES_bl ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(bl.mod.full)
step(bl.mod.full)

bl.mod.null <- glm(SES_bl ~ 1, family = gaussian, data = SES)
summary(bl.mod.null)
qqnorm(resid(bl.mod.null))
qqline(resid(bl.mod.null))
plot(simulateResiduals(bl.mod.null))
densityPlot(rstudent(bl.mod.null)) # check density estimate of the distribution of residuals
outlierTest(bl.mod.null)
influenceIndexPlot(bl.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(bl.mod.full, bl.mod.null, test = "F")
AICctab(bl.mod.full, bl.mod.null)



## lecty---- 
##lec_0- Kleptoparasitic
hist(SES$SES_lec_0)
plot(SES$SES_lec_0, ylim = c(-3.5, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_lec_0 ~ SES$trmt, ylim = c(-3.5, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
lec_0.farm <- wilcox.test(farm$SES_lec_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_0.farm
lec_0.control <- wilcox.test(control$SES_lec_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_0.control
lec_0.t1 <- wilcox.test(T1$SES_lec_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_0.t1
lec_0.t8 <- wilcox.test(T8$SES_lec_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_0.t8


## compare among treatments and landscape variables
dotchart(SES$SES_lec_0, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_lec_0 ~ trmt))
with(SES, ad.test(SES_lec_0))

lec_0.mod.full <- glm(SES_lec_0 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(lec_0.mod.full)
step(lec_0.mod.full)

lec_0.mod.null <- glm(SES_lec_0 ~ 1, family = gaussian, data = SES)
summary(lec_0.mod.null)
qqnorm(resid(lec_0.mod.null))
qqline(resid(lec_0.mod.null))
plot(simulateResiduals(lec_0.mod.null))
densityPlot(rstudent(lec_0.mod.null)) # check density estimate of the distribution of residuals
outlierTest(lec_0.mod.null)
influenceIndexPlot(lec_0.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(lec_0.mod.full, lec_0.mod.null, test = "F")
AICctab(lec_0.mod.full, lec_0.mod.null)



## lec_1 - Generalist
hist(SES$SES_lec_1)
SES$SES_lec_1
plot(SES$SES_lec_1)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_lec_1 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
lec_1.farm <- wilcox.test(farm$SES_lec_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_1.farm
lec_1.control <- wilcox.test(control$SES_lec_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_1.control
lec_1.t1 <- wilcox.test(T1$SES_lec_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_1.t1
lec_1.t8 <- wilcox.test(T8$SES_lec_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_1.t8

## compare among treatments and landscape variables
dotchart(SES$SES_lec_1, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_lec_1 ~ trmt))
with(SES, ad.test(SES_lec_1))

lec_1.mod.full <- glm(SES_lec_1 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(lec_1.mod.full)
step(lec_1.mod.full)

lec_1.mod.red <- glm(SES_lec_1 ~ trmt, family = gaussian, data = SES)
summary(lec_1.mod.red)
qqnorm(resid(lec_1.mod.red))
qqline(resid(lec_1.mod.red))
plot(simulateResiduals(lec_1.mod.red))
densityPlot(rstudent(lec_1.mod.red)) # check density estimate of the distribution of residuals
outlierTest(lec_1.mod.red)
influenceIndexPlot(lec_1.mod.red, vars = c("Cook"), id = list(n = 3))

# significant outlier, let's remove it and see if the model fits better
lec_1.mod.red2 <- update(lec_1.mod.red, subset = -c(1))
summary(lec_1.mod.red2)
compareCoefs(lec_1.mod.red, lec_1.mod.red2) # compares estimated coefficients and their standard errors

Anova(lec_1.mod.red)
emmeans(lec_1.mod.red, pairwise ~ trmt)

effect_plot(lec_1.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'Largest Patch Index (LPI)', y.label = 'Standardized Effect Sizes (SES)')

lec_1.mod.null <- glm(SES_lec_1 ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(lec_1.mod.full, lec_1.mod.red, lec_1.mod.null, test = "F")
AICctab(lec_1.mod.full, lec_1.mod.red, lec_1.mod.null)



## lec_2 - Specialist
hist(SES$SES_lec_2)
SES$SES_lec_2
plot(SES$SES_lec_2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_lec_2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
lec_2.farm <- wilcox.test(farm$SES_lec_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_2.farm
lec_2.control <- wilcox.test(control$SES_lec_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_2.control
lec_2.t1 <- wilcox.test(T1$SES_lec_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_2.t1
lec_2.t8 <- wilcox.test(T8$SES_lec_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_2.t8


## compare among treatments and landscape variables
dotchart(SES$SES_lec_2, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_lec_2 ~ trmt))
with(SES, ad.test(SES_lec_2))

lec_2.mod.full <- glm(SES_lec_2 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(lec_2.mod.full)
step(lec_2.mod.full)

lec_2.mod.red <- glm(SES_lec_2 ~ trmt + lpi, family = gaussian, data = SES)
summary(lec_2.mod.red)
qqnorm(resid(lec_2.mod.red))
qqline(resid(lec_2.mod.red))
plot(simulateResiduals(lec_2.mod.red))
densityPlot(rstudent(lec_2.mod.red)) # check density estimate of the distribution of residuals
outlierTest(lec_2.mod.red)
influenceIndexPlot(lec_2.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(lec_2.mod.red)
emmeans(lec_2.mod.red, pairwise ~ trmt)
plot_summs(lec_2.mod.red, scale = TRUE)

effect_plot(lec_2.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(lec_2.mod.red, pred = lpi, interval = TRUE, partial.residuals = TRUE, x.label = 'Largest Patch Index (LPI)', y.label = 'Standardized Effect Sizes (SES)')

lec_2.mod.null <- glm(SES_lec_2 ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(lec_2.mod.full, lec_2.mod.red, lec_2.mod.null, test = "F")
AICctab(lec_2.mod.full, lec_2.mod.red, lec_2.mod.null)


##Nesting----
## nest_1 - Soil
hist(SES$SES_nest_1)
plot(SES$SES_nest_1)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest_1 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
nest_1.farm <- wilcox.test(farm$SES_nest_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_1.farm
nest_1.control <- wilcox.test(control$SES_nest_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_1.control
nest_1.t1 <- wilcox.test(T1$SES_nest_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_1.t1
nest_1.t8 <- wilcox.test(T8$SES_nest_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_1.t8


## compare among treatments and landscape variables
dotchart(SES$SES_nest_1, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_nest_1 ~ trmt))
with(SES, ad.test(SES_nest_1))

nest_1.mod.full <- glm(SES_nest_1 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(nest_1.mod.full)
step(nest_1.mod.full)

nest_1.mod.red <- glm(SES_nest_1 ~ trmt + pland + enn, family = gaussian, data = SES)
summary(nest_1.mod.red)
#the below code is for seeing if variables are correlated
summ(nest_1.mod.red, scale = TRUE, confint = TRUE, vifs = TRUE) #VIF should be <3
qqnorm(resid(nest_1.mod.red))
qqline(resid(nest_1.mod.red))
plot(simulateResiduals(nest_1.mod.red))
densityPlot(rstudent(nest_1.mod.red)) # check density estimate of the distribution of residuals
outlierTest(nest_1.mod.red)
influenceIndexPlot(nest_1.mod.red, vars = c("Cook"), id = list(n = 3))

# significant outlier, let's remove it and see if it improves model fit
nest_1.mod.red2 <- update(nest_1.mod.red, subset = -c(3))
summary(nest_1.mod.red2)
compareCoefs(nest_1.mod.red, nest_1.mod.red2) # compares estimated coefficients and their standard errors

Anova(nest_1.mod.red)
emmeans(nest_1.mod.red, pairwise ~ trmt)

effect_plot(nest_1.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(nest_1.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = 'Percentage Greenspace', y.label = 'Standardized Effect Sizes (SES)')

nest_1.mod.null <- glm(SES_nest_1 ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(nest_1.mod.full, nest_1.mod.red, nest_1.mod.null)
AICctab(nest_1.mod.full, nest_1.mod.red, nest_1.mod.null)



## nest_2 - Cavity
hist(SES$SES_nest_2)
plot(SES$SES_nest_2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest_2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
nest_2.farm <- wilcox.test(farm$SES_nest_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_2.farm
nest_2.control <- wilcox.test(control$SES_nest_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_2.control
nest_2.t1 <- wilcox.test(T1$SES_nest_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_2.t1
nest_2.t8 <- wilcox.test(T8$SES_nest_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_2.t8


## compare among treatments and landscape variables
dotchart(SES$SES_nest_2, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_nest_2 ~ trmt))
with(SES, ad.test(SES_nest_2))

nest_2.mod.full <- glm(SES_nest_2 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(nest_2.mod.full)
step(nest_2.mod.full)

nest_2.mod.red <- glm(SES_nest_2 ~ trmt, family = gaussian, data = SES)
summary(nest_2.mod.red)
qqnorm(resid(nest_2.mod.red))
qqline(resid(nest_2.mod.red))
plot(simulateResiduals(nest_2.mod.red))
densityPlot(rstudent(nest_2.mod.red)) # check density estimate of the distribution of residuals
outlierTest(nest_2.mod.red)
influenceIndexPlot(nest_2.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(nest_2.mod.red)
emmeans(nest_2.mod.red, pairwise ~ trmt)

effect_plot(nest_2.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')

nest_2.mod.null <- glm(SES_nest_2 ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(nest_2.mod.full, nest_2.mod.red, nest_2.mod.null)
AICctab(nest_2.mod.full, nest_2.mod.red, nest_2.mod.null)



## nest_3 - Hive
hist(SES$SES_nest_3)
plot(SES$SES_nest_3)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest_3 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
nest_3.farm <- wilcox.test(farm$SES_nest_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_3.farm
nest_3.control <- wilcox.test(control$SES_nest_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_3.control
nest_3.t1 <- wilcox.test(T1$SES_nest_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_3.t1
nest_3.t8 <- wilcox.test(T8$SES_nest_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_3.t8


## compare among treatments and landscape variables
dotchart(SES$SES_nest_3, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_nest_3 ~ trmt))
with(SES, ad.test(SES_nest_3))

nest_3.mod.full <- glm(SES_nest_3 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(nest_3.mod.full)
step(nest_3.mod.full)

nest_3.mod.red <- glm(SES_nest_3 ~ trmt, family = gaussian, data = SES)
summary(nest_3.mod.red)
qqnorm(resid(nest_3.mod.red))
qqline(resid(nest_3.mod.red))
plot(simulateResiduals(nest_3.mod.red))
densityPlot(rstudent(nest_3.mod.red)) # check density estimate of the distribution of residuals
outlierTest(nest_3.mod.red)
influenceIndexPlot(nest_3.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(nest_3.mod.red)
emmeans(nest_3.mod.red, pairwise ~ trmt)

effect_plot(nest_3.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')

nest_3.mod.null <- glm(SES_nest_3 ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(nest_3.mod.full, nest_3.mod.red, nest_3.mod.null)
AICctab(nest_3.mod.full, nest_3.mod.red, nest_3.mod.null)



## nest_4 - Pithy Stems
hist(SES$SES_nest_4)
plot(SES$SES_nest_4)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest_4 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
nest_4.farm <- wilcox.test(farm$SES_nest_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_4.farm
nest_4.control <- wilcox.test(control$SES_nest_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_4.control
nest_4.t1 <- wilcox.test(T1$SES_nest_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_4.t1
nest_4.t8 <- wilcox.test(T8$SES_nest_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_4.t8


## compare among treatments and landscape variables
dotchart(SES$SES_nest_4, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_nest_4 ~ trmt))
with(SES, ad.test(SES_nest_4))

nest_4.mod.full <- glm(SES_nest_4 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(nest_4.mod.full)
step(nest_4.mod.full)

summary(nest_4.mod.full)
summ(nest_4.mod.full, scale = TRUE, confint = TRUE, vifs = TRUE)
qqnorm(resid(nest_4.mod.full))
qqline(resid(nest_4.mod.full))
plot(simulateResiduals(nest_4.mod.full))
densityPlot(rstudent(nest_4.mod.full)) # check density estimate of the distribution of residuals
outlierTest(nest_4.mod.full)
influenceIndexPlot(nest_4.mod.full, vars = c("Cook"), id = list(n = 3))

Anova(nest_4.mod.full)

effect_plot(nest_4.mod.full, pred = lpi, interval = TRUE, partial.residuals = TRUE, x.label = 'Largest Patch Index', y.label = 'Standardized Effect Sizes (SES)')

nest_4.mod.null <- glm(SES_nest_4 ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(nest_4.mod.full, nest_4.mod.null)
AICctab(nest_4.mod.full, nest_4.mod.null)



## nest_5 - Wood
hist(SES$SES_nest_5)
plot(SES$SES_nest_5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest_5 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
nest_5.farm <- wilcox.test(farm$SES_nest_5, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_5.farm
nest_5.control <- wilcox.test(control$SES_nest_5, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_5.control
nest_5.t1 <- wilcox.test(T1$SES_nest_5, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_5.t1
nest_5.t8 <- wilcox.test(T8$SES_nest_5, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_5.t8


## compare among treatments and landscape variables
dotchart(SES$SES_nest_5, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_nest_5 ~ trmt))
with(SES, ad.test(SES_nest_5))

nest_5.mod.full <- glm(SES_nest_5 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(nest_5.mod.full)
step(nest_5.mod.full)

nest_5.mod.red <- glm(SES_nest_5 ~ trmt + lpi, family = gaussian, data = SES)
summary(nest_5.mod.red)
qqnorm(resid(nest_5.mod.red))
qqline(resid(nest_5.mod.red))
plot(simulateResiduals(nest_5.mod.red))
densityPlot(rstudent(nest_5.mod.red)) # check density estimate of the distribution of residuals
outlierTest(nest_5.mod.red)
influenceIndexPlot(nest_5.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(nest_5.mod.red)
emmeans(nest_5.mod.red, pairwise ~ trmt)

effect_plot(nest_5.mod.full, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'Largest Patch Index', y.label = 'Standardized Effect Sizes (SES)')

nest_5.mod.null <- glm(SES_nest_5 ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(nest_5.mod.full, nest_5.mod.red, nest_5.mod.null)
AICctab(nest_5.mod.full, nest_5.mod.red, nest_5.mod.null)


##Sociality----
## soc_1 - Subsocial
hist(SES$SES_soc_1)
plot(SES$SES_soc_1)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_soc_1 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
soc_1.farm <- wilcox.test(farm$SES_soc_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_1.farm
soc_1.control <- wilcox.test(control$SES_soc_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_1.control
soc_1.t1 <- wilcox.test(T1$SES_soc_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_1.t1
soc_1.t8 <- wilcox.test(T8$SES_soc_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_1.t8


## compare among treatments and landscape variables
dotchart(SES$SES_soc_1, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_soc_1 ~ trmt))
with(SES, ad.test(SES_soc_1))

soc_1.mod.full <- glm(SES_soc_1 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(soc_1.mod.full)
step(soc_1.mod.full)

soc_1.mod.red <- glm(SES_soc_1 ~ pland + enn, family = gaussian, data = SES)
summary(soc_1.mod.red)
summ(soc_1.mod.red, scale = TRUE, confint = TRUE, vifs = TRUE)
qqnorm(resid(soc_1.mod.red))
qqline(resid(soc_1.mod.red))
plot(simulateResiduals(soc_1.mod.red))
densityPlot(rstudent(soc_1.mod.red)) # check density estimate of the distribution of residuals
outlierTest(soc_1.mod.red)
influenceIndexPlot(soc_1.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(soc_1.mod.red)

effect_plot(soc_1.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = 'Percentge Greenspace', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(soc_1.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Greenspace Isolation (ENN)', y.label = 'Standardized Effect Sizes (SES)')

soc_1.mod.null <- glm(SES_soc_1 ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(soc_1.mod.full, soc_1.mod.red, soc_1.mod.null)
AICctab(soc_1.mod.full, soc_1.mod.red, soc_1.mod.null)



## soc_2 - Solitary
hist(SES$SES_soc_2)
plot(SES$SES_soc_2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_soc_2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
soc_2.farm <- wilcox.test(farm$SES_soc_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_2.farm
soc_2.control <- wilcox.test(control$SES_soc_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_2.control
soc_2.t1 <- wilcox.test(T1$SES_soc_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_2.t1
soc_2.t8 <- wilcox.test(T8$SES_soc_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_2.t8


## compare among treatments and landscape variables
dotchart(SES$SES_soc_2, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_soc_2 ~ trmt))
with(SES, ad.test(SES_soc_2))

soc_2.mod.full <- glm(SES_soc_2 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(soc_2.mod.full)
step(soc_2.mod.full)

soc_2.mod.red <- glm(SES_soc_2 ~ lpi, family = gaussian, data = SES)
summary(soc_2.mod.red)
qqnorm(resid(soc_2.mod.red))
qqline(resid(soc_2.mod.red))
plot(simulateResiduals(soc_2.mod.red))
densityPlot(rstudent(soc_2.mod.red)) # check density estimate of the distribution of residuals
outlierTest(soc_2.mod.red)
influenceIndexPlot(soc_2.mod.red, vars = c("Cook"), id = list(n = 3))

# significant outlier, let's remove it and see if the model fits better
soc_2.mod.red2 <- update(soc_2.mod.red, subset = -c(13))
summary(soc_2.mod.red2)
compareCoefs(soc_2.mod.red, soc_2.mod.red2) # compares estimated coefficients and their standard errors

Anova(soc_2.mod.red)

effect_plot(soc_2.mod.red, pred = lpi, interval = TRUE, partial.residuals = TRUE, x.label = 'Largest Patch Index (LPI)', y.label = 'Standardized Effect Sizes (SES)')

soc_2.mod.null <- glm(SES_soc_2 ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(soc_2.mod.full, soc_2.mod.red, soc_2.mod.null)
AICctab(soc_2.mod.full, soc_2.mod.red, soc_2.mod.null)



## soc_3 - Eusocial
hist(SES$SES_soc_3)
plot(SES$SES_soc_3)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_soc_3 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
soc_3.farm <- wilcox.test(farm$SES_soc_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_3.farm
soc_3.control <- wilcox.test(control$SES_soc_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_3.control
soc_3.t1 <- wilcox.test(T1$SES_soc_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_3.t1
soc_3.t8 <- wilcox.test(T8$SES_soc_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_3.t8


## compare among treatments and landscape variables
dotchart(SES$SES_soc_3, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_soc_3 ~ trmt))
with(SES, ad.test(SES_soc_3))

soc_3.mod.full <- glm(SES_soc_3 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(soc_3.mod.full)
step(soc_3.mod.full)


soc_3.mod.null <- glm(SES_soc_3 ~ 1, family = gaussian, data = SES)
summary(soc_3.mod.null)
qqnorm(resid(soc_3.mod.null))
qqline(resid(soc_3.mod.null))
plot(simulateResiduals(soc_3.mod.null))
densityPlot(rstudent(soc_3.mod.null)) # check density estimate of the distribution of residuals
outlierTest(soc_3.mod.null)
influenceIndexPlot(soc_3.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(soc_3.mod.full, soc_3.mod.null)
AICctab(soc_3.mod.full, soc_3.mod.null)



## soc_4 - Parasitic
hist(SES$SES_soc_4)
plot(SES$SES_soc_4, ylim = c(-3.5, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_soc_4 ~ SES$trmt, ylim = c(-3.5, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
soc_4.farm <- wilcox.test(farm$SES_soc_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_4.farm
soc_4.control <- wilcox.test(control$SES_soc_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_4.control
soc_4.t1 <- wilcox.test(T1$SES_soc_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_4.t1
soc_4.t8 <- wilcox.test(T8$SES_soc_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_4.t8


## compare among treatments and landscape variables
dotchart(SES$SES_soc_4, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_soc_4 ~ trmt))
with(SES, ad.test(SES_soc_4))

soc_4.mod.full <- glm(SES_soc_4 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(soc_4.mod.full)
step(soc_4.mod.full)

soc_4.mod.null <- glm(SES_soc_4 ~ 1, family = gaussian, data = SES)
summary(soc_4.mod.null)
qqnorm(resid(soc_4.mod.null))
qqline(resid(soc_4.mod.null))
plot(simulateResiduals(soc_4.mod.null))
densityPlot(rstudent(soc_4.mod.null)) # check density estimate of the distribution of residuals
outlierTest(soc_4.mod.null)
influenceIndexPlot(soc_4.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(soc_4.mod.full, soc_4.mod.null)
AICctab(soc_4.mod.full, soc_4.mod.null)


##Origin----
## ori_0 - Native
hist(SES$SES_ori_0)
plot(SES$SES_ori_0)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_ori_0 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
ori_0.farm <- wilcox.test(farm$SES_ori_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_0.farm
ori_0.control <- wilcox.test(control$SES_ori_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_0.control
ori_0.t1 <- wilcox.test(T1$SES_ori_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_0.t1
ori_0.t8 <- wilcox.test(T8$SES_ori_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_0.t8


## compare among treatments and landscape variables
dotchart(SES$SES_ori_0, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_ori_0 ~ trmt))
with(SES, ad.test(SES_ori_0))

ori_0.mod.full <- glm(SES_ori_0 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(ori_0.mod.full)
step(ori_0.mod.full)

ori_0.mod.red <- glm(SES_ori_0 ~ trmt, family = gaussian, data = SES)
summary(ori_0.mod.red)
qqnorm(resid(ori_0.mod.red))
qqline(resid(ori_0.mod.red))
plot(simulateResiduals(ori_0.mod.red))
densityPlot(rstudent(ori_0.mod.red)) # check density estimate of the distribution of residuals
outlierTest(ori_0.mod.red)
influenceIndexPlot(ori_0.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(ori_0.mod.red)

effect_plot(ori_0.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'Largest Patch Index (LPI)', y.label = 'Standardized Effect Sizes (SES)')

ori_0.mod.null <- glm(SES_ori_0 ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(ori_0.mod.full, ori_0.mod.red, ori_0.mod.null)
AICctab(ori_0.mod.full, ori_0.mod.red, ori_0.mod.null)



## ori_2 - Exotic
hist(SES$SES_ori_1)
plot(SES$SES_ori_1)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_ori_1 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
ori_1.farm <- wilcox.test(farm$SES_ori_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_1.farm
ori_1.control <- wilcox.test(control$SES_ori_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_1.control
ori_1.t1 <- wilcox.test(T1$SES_ori_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_1.t1
ori_1.t8 <- wilcox.test(T8$SES_ori_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_1.t8


## compare among treatments and landscape variables
dotchart(SES$SES_ori_1, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_ori_1 ~ trmt))
with(SES, ad.test(SES_ori_1))

ori_1.mod.full <- glm(SES_ori_1 ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(ori_1.mod.full)
step(ori_1.mod.full)

ori_1.mod.red <- glm(SES_ori_1 ~ trmt, family = gaussian, data = SES)
summary(ori_1.mod.red)
qqnorm(resid(ori_1.mod.red))
qqline(resid(ori_1.mod.red))
plot(simulateResiduals(ori_1.mod.red))
densityPlot(rstudent(ori_1.mod.red)) # check density estimate of the distribution of residuals
outlierTest(ori_1.mod.red)
influenceIndexPlot(ori_1.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(ori_1.mod.red)

effect_plot(ori_1.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'Largest Patch Index (LPI)', y.label = 'Standardized Effect Sizes (SES)')

ori_1.mod.null <- glm(SES_ori_1 ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(ori_1.mod.full, ori_1.mod.red, ori_1.mod.null)
AICctab(ori_1.mod.full, ori_1.mod.red, ori_1.mod.null)


##taxonomic diversity----
## taxonomic diveristy - beta sor
hist(SES$SES_bsor)
plot(SES$SES_bsor, ylim = c(-0.5, 8))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsor ~ SES$trmt, ylim = c(-0.5, 8))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
bsor.farm <- wilcox.test(farm$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsor.farm
bsor.control <- wilcox.test(control$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsor.control
bsor.t1 <- wilcox.test(T1$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsor.t1
bsor.t8 <- wilcox.test(T8$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsor.t8


## compare among treatments and landscape variables
dotchart(SES$SES_bsor, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_bsor ~ trmt))
with(SES, ad.test(SES_bsor))

bsor.mod.full <- glm(SES_bsor ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(bsor.mod.full)
step(bsor.mod.full)

bsor.mod.red <- glm(SES_bsor ~ enn, family = gaussian, data = SES)
summary(bsor.mod.red)
qqnorm(resid(bsor.mod.red))
qqline(resid(bsor.mod.red))
plot(simulateResiduals(bsor.mod.red))
densityPlot(rstudent(bsor.mod.red)) # check density estimate of the distribution of residuals
outlierTest(bsor.mod.red)
influenceIndexPlot(bsor.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(bsor.mod.red)

effect_plot(bsor.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Greenspace Isolation (ENN)', y.label = 'Standardized Effect Sizes (SES)')

bsor.mod.null <- glm(SES_bsor ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(bsor.mod.full, bsor.mod.red, bsor.mod.null)
AICctab(bsor.mod.full, bsor.mod.red, bsor.mod.null)



## taxonomic diveristy - beta sim
hist(SES$SES_bsim)
plot(SES$SES_bsim)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
bsim.farm <- wilcox.test(farm$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsim.farm
bsim.control <- wilcox.test(control$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsim.control
bsim.t1 <- wilcox.test(T1$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsim.t1
bsim.t8 <- wilcox.test(T8$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsim.t8


## compare among treatments and landscape variables
dotchart(SES$SES_bsim, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_bsim ~ trmt))
with(SES, ad.test(SES_bsim))

bsim.mod.full <- glm(SES_bsim ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(bsim.mod.full)
step(bsim.mod.full)

bsim.mod.red <- glm(SES_bsim ~ enn, family = gaussian, data = SES)
summary(bsim.mod.red)
qqnorm(resid(bsim.mod.red))
qqline(resid(bsim.mod.red))
plot(simulateResiduals(bsim.mod.red))
densityPlot(rstudent(bsim.mod.red)) # check density estimate of the distribution of residuals
outlierTest(bsim.mod.red)
influenceIndexPlot(bsim.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(bsim.mod.red)

effect_plot(bsim.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Greenspace Isolation (ENN)', y.label = 'Standardized Effect Sizes (SES)')

bsim.mod.null <- glm(SES_bsim ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(bsim.mod.full, bsim.mod.red, bsim.mod.null)
AICctab(bsim.mod.full, bsim.mod.red, bsim.mod.null)



## taxonomic diveristy - beta sne
hist(SES$SES_bsne)
plot(SES$SES_bsne)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsne ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
bsne.farm <- wilcox.test(farm$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsne.farm
bsne.control <- wilcox.test(control$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsne.control
bsne.t1 <- wilcox.test(T1$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsne.t1
bsne.t8 <- wilcox.test(T8$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsne.t8


## compare among treatments and landscape variables
dotchart(SES$SES_bsne, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_bsne ~ trmt))
with(SES, ad.test(SES_bsne))

bsne.mod.full <- glm(SES_bsne ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(bsne.mod.full)
step(bsne.mod.full)

bsne.mod.null <- glm(SES_bsne ~ 1, family = gaussian, data = SES)
summary(bsne.mod.null)
qqnorm(resid(bsne.mod.null))
qqline(resid(bsne.mod.null))
plot(simulateResiduals(bsne.mod.null))
densityPlot(rstudent(bsne.mod.null)) # check density estimate of the distribution of residuals
outlierTest(bsne.mod.null)
influenceIndexPlot(bsne.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(bsne.mod.full, bsne.mod.null)
AICctab(bsne.mod.full, bsne.mod.null)



## functional alpha diversity----
hist(SES$SES_falpha)
plot(SES$SES_falpha)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_falpha ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
falpha.farm <- wilcox.test(farm$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
falpha.farm
falpha.control <- wilcox.test(control$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
falpha.control
falpha.t1 <- wilcox.test(T1$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
falpha.t1
falpha.t8 <- wilcox.test(T8$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
falpha.t8


## compare among treatments and landscape variables
dotchart(SES$SES_falpha, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_falpha ~ trmt))
with(SES, ad.test(SES_falpha))

falpha.mod.full <- glm(SES_falpha ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(falpha.mod.full)
step(falpha.mod.full)

falpha.mod.red <- glm(SES_falpha ~ trmt, family = gaussian, data = SES)
summary(falpha.mod.red)
qqnorm(resid(falpha.mod.red))
qqline(resid(falpha.mod.red))
plot(simulateResiduals(falpha.mod.red))
densityPlot(rstudent(falpha.mod.red)) # check density estimate of the distribution of residuals
outlierTest(falpha.mod.red)
influenceIndexPlot(falpha.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(falpha.mod.red)
emmeans(falpha.mod.red, pairwise ~ trmt)

effect_plot(falpha.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')

falpha.mod.null <- glm(SES_falpha ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(falpha.mod.full, falpha.mod.red, falpha.mod.null)
AICtab(falpha.mod.full, falpha.mod.red, falpha.mod.null)


##Functional beta div----
## functional beta diversity - beta sor
hist(SES$SES_fbsor)
plot(SES$SES_fbsor, ylim = c(-0.5, 6.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsor ~ SES$trmt, ylim = c(-0.5, 6.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
fbsor.farm <- wilcox.test(farm$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsor.farm
fbsor.control <- wilcox.test(control$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsor.control
fbsor.t1 <- wilcox.test(T1$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsor.t1
fbsor.t8 <- wilcox.test(T8$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsor.t8


## compare among treatments and landscape variables
dotchart(SES$SES_fbsor, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_fbsor ~ trmt))
with(SES, ad.test(SES_fbsor))

fbsor.mod.full <- glm(SES_fbsor ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(fbsor.mod.full)
step(fbsor.mod.full)

fbsor.mod.red <- glm(SES_fbsor ~ trmt, family = gaussian, data = SES)
summary(fbsor.mod.red)
qqnorm(resid(fbsor.mod.red))
qqline(resid(fbsor.mod.red))
plot(simulateResiduals(fbsor.mod.red))
densityPlot(rstudent(fbsor.mod.red)) # check density estimate of the distribution of residuals
outlierTest(fbsor.mod.red)
influenceIndexPlot(fbsor.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(fbsor.mod.red)
emmeans(fbsor.mod.red, pairwise ~ trmt)

effect_plot(fbsor.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')

fbsor.mod.null <- glm(SES_fbsor ~ 1, family = gaussian, data = SES)

# model comparison techniques
anova(fbsor.mod.full, fbsor.mod.red, fbsor.mod.null)
AICctab(fbsor.mod.full, fbsor.mod.red, fbsor.mod.null)



## functional beta diversity - beta sim
hist(SES$SES_fbsim)
plot(SES$SES_fbsim)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
fbsim.farm <- wilcox.test(farm$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsim.farm
fbsim.control <- wilcox.test(control$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsim.control
fbsim.t1 <- wilcox.test(T1$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsim.t1
fbsim.t8 <- wilcox.test(T8$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsim.t8


## compare among treatments and landscape variables
dotchart(SES$SES_fbsim, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_fbsim ~ trmt))
with(SES, ad.test(SES_fbsim))
 
fbsim.mod.full <- glm(SES_fbsim ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(fbsim.mod.full)
step(fbsim.mod.full)

fbsim.mod.null <- glm(SES_fbsim ~ 1, family = gaussian, data = SES)
summary(fbsim.mod.null)
qqnorm(resid(fbsim.mod.null))
qqline(resid(fbsim.mod.null))
plot(simulateResiduals(fbsim.mod.null))
densityPlot(rstudent(fbsim.mod.null)) # check density estimate of the distribution of residuals
outlierTest(fbsim.mod.null)
influenceIndexPlot(fbsim.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(fbsim.mod.full, fbsim.mod.null)
AICctab(fbsim.mod.full, fbsim.mod.null)



## functional beta diversity - beta sne
hist(SES$SES_fbsne)
plot(SES$SES_fbsne, ylim = c(-0.5, 6))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsne ~ SES$trmt, ylim = c(-0.5, 6))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
fbsne.farm <- wilcox.test(farm$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsne.farm
fbsne.control <- wilcox.test(control$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsne.control
fbsne.t1 <- wilcox.test(T1$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsne.t1
fbsne.t8 <- wilcox.test(T8$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsne.t8


## compare among treatments and landscape variables
dotchart(SES$SES_fbsne, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_fbsne ~ trmt))
with(SES, ad.test(SES_fbsne))

fbsne.mod.full <- glm(SES_fbsne ~ trmt + pland + lpi + enn, family = gaussian, data = SES)
summary(fbsne.mod.full)
step(fbsne.mod.full)

fbsne.mod.null <- glm(SES_fbsne ~ 1, family = gaussian, data = SES)
summary(fbsne.mod.null)
qqnorm(resid(fbsne.mod.null))
qqline(resid(fbsne.mod.null))
plot(simulateResiduals(fbsne.mod.null))
densityPlot(rstudent(fbsne.mod.null)) # check density estimate of the distribution of residuals
outlierTest(fbsne.mod.null)
influenceIndexPlot(fbsne.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(fbsne.mod.full, fbsne.mod.null)
AICctab(fbsne.mod.full, fbsne.mod.null)







