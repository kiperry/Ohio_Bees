###################################################################################
#
# Ohio bee data- FS
#
# Step 2: Urban Cleveland Pool to Local Greenspaces
# Comparing among Local habitats in Frances Sivakoff's data
# 
#CWM and Functional Diversity - Null Models
# GLMMs
#
# 
#CA Shepard and KI Perry: 2 October 2023
###################################################################################

#Refer back to Final_local_analysis to get to the beginning of this document
#Everything before Line 83
t <- read.csv("traits_urbanpool.csv", row.names=1)
a <- read.csv("urbanpool.FS.csv", row.names=1)

aO<-a
a<-a[1:134]

str(t)
names(t)
colnames(t) <- c("bl", "lec", "nest", "soc", "ori")
names(t)



plot(t)
cor(t, method = c("pearson"), use = "complete.obs")
#.6 is a bit high, but not going to remove sociality or lecty
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

#just testing file locations before I continue
write.csv(t, file = "Urban to Local_Nulls_FS/test.csv")
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
wt <- c(0.39, 0.18, 0.11, 0.12, 0.20)

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
  nsp <- randomizeMatrix(samp = a, null.model = "richness")#richness means it will accept the 0s
  
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
write.csv(nbl, file = "Urban to Local_Nulls_FS/nbl_FS.csv")
write.csv(nlec_0, file = "Urban to Local_Nulls_FS/nlec_0_FS.csv")
write.csv(nlec_1, file = "Urban to Local_Nulls_FS/nlec_1_FS.csv")
write.csv(nlec_2, file = "Urban to Local_Nulls_FS/nlec_2_FS.csv")
write.csv(nnest_1, file = "Urban to Local_Nulls_FS/nnest_1_FS.csv")
write.csv(nnest_2, file = "Urban to Local_Nulls_FS/nnest_2_FS.csv")
write.csv(nnest_3, file = "Urban to Local_Nulls_FS/nnest_3_FS.csv")
write.csv(nnest_4, file = "Urban to Local_Nulls_FS/nnest_4_FS.csv")
write.csv(nnest_5, file = "Urban to Local_Nulls_FS/nnest_5_FS.csv")
write.csv(nsoc_1, file = "Urban to Local_Nulls_FS/nsoc_1_FS.csv")
write.csv(nsoc_2, file = "Urban to Local_Nulls_FS/nsoc_2_FS.csv")
write.csv(nsoc_3, file = "Urban to Local_Nulls_FS/nsoc_3_FS.csv")
write.csv(nsoc_4, file = "Urban to Local_Nulls_FS/nsoc_4_FS.csv")
write.csv(nori_0, file = "Urban to Local_Nulls_FS/nori_0_FS.csv")
write.csv(nori_1, file = "Urban to Local_Nulls_FS/nori_1_FS.csv")

write.csv(nbsim, file = "Urban to Local_Nulls_FS/tbeta_sim_FS.csv")
write.csv(nbsne, file = "Urban to Local_Nulls_FS/tbeta_sne_FS.csv")
write.csv(nbsor, file = "Urban to Local_Nulls_FS/tbeta_sor_FS.csv")

write.csv(nfalpha, file = "Urban to Local_Nulls_FS/falpha_FS.csv")

write.csv(nfsim, file = "Urban to Local_Nulls_FS/fbeta_sim_FS.csv")
write.csv(nfsne, file = "Urban to Local_Nulls_FS/fbeta_sne_FS.csv")
write.csv(nfsor, file = "Urban to Local_Nulls_FS/fbeta_sor_FS.csv")

# load the output matrices

nbl <- read.csv("Urban to Local_Nulls_FS/nbl_FS.csv", row.names=1)
nlec_0 <- read.csv("Urban to Local_Nulls_FS/nlec_0_FS.csv", row.names=1)
nlec_1 <- read.csv("Urban to Local_Nulls_FS/nlec_1_FS.csv", row.names=1)
nlec_2 <- read.csv("Urban to Local_Nulls_FS/nlec_2_FS.csv", row.names=1)
nnest_1 <- read.csv("Urban to Local_Nulls_FS/nnest_1_FS.csv", row.names=1)
nnest_2 <- read.csv("Urban to Local_Nulls_FS/nnest_2_FS.csv", row.names=1)
nnest_3 <- read.csv("Urban to Local_Nulls_FS/nnest_3_FS.csv", row.names=1)
nnest_4 <- read.csv("Urban to Local_Nulls_FS/nnest_4_FS.csv", row.names=1)
nnest_5 <- read.csv("Urban to Local_Nulls_FS/nnest_5_FS.csv", row.names=1)
nsoc_1 <- read.csv("Urban to Local_Nulls_FS/nsoc_1_FS.csv", row.names=1)
nsoc_2 <- read.csv("Urban to Local_Nulls_FS/nsoc_2_FS.csv", row.names=1)
nsoc_3 <- read.csv("Urban to Local_Nulls_FS/nsoc_3_FS.csv", row.names=1)
nsoc_4 <- read.csv("Urban to Local_Nulls_FS/nsoc_4_FS.csv", row.names=1)
nori_0 <- read.csv("Urban to Local_Nulls_FS/nori_0_FS.csv", row.names=1)
nori_1 <- read.csv("Urban to Local_Nulls_FS/nori_1_FS.csv", row.names=1)

nfalpha <- read.csv("Urban to Local_Nulls_FS/falpha_FS.csv", row.names=1)

nbsim <- read.csv("Urban to Local_Nulls_FS/tbeta_sim_FS.csv", row.names=1)
nbsne <- read.csv("Urban to Local_Nulls_FS/tbeta_sne_FS.csv", row.names=1)
nbsor <- read.csv("Urban to Local_Nulls_FS/tbeta_sim_FS.csv", row.names=1)

nfsim <- read.csv("Urban to Local_Nulls_FS/fbeta_sim_FS.csv", row.names=1)
nfsne <- read.csv("Urban to Local_Nulls_FS/fbeta_sne_FS.csv", row.names=1)
nfsor <- read.csv("Urban to Local_Nulls_FS/fbeta_sor_FS.csv", row.names=1)

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

write.csv(SES.all, file = "Urban to Local_Nulls_FS/SES_Local_FS.csv")
#import the SES data
SES <- read.csv("Urban to Local_Nulls_FS/SES_Local_FS.csv", row.names = 1)
SES

#Including treatment

SES$trmt <- aO$trmt
str(SES)

#Adding the random site variable
SES$site <- c("47thStGarden", "BluePike", "Buckeye", "Esperanza", "Fairfax","LonnieBurten","Midtown","Slavicvillage","C10","C11","C13","C5","C6","C69","C7","C9")


## pull out data for each treatment
farm <- SES[which(SES$trmt == "Farm"),]
str(farm)

control <- SES[which(SES$trmt == "Control"),]
str(control)

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


##taxonomic diversity----
### taxonomic diveristy - beta sor----
hist(SES$SES_bsor)
plot(SES$SES_bsor, ylim = c(-30, 8))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsor ~ SES$trmt, ylim = c(-19, 1))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

# compare to null expectations by treatment
bsor.farm <- wilcox.test(farm$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsor.farm
bsor.control <- wilcox.test(control$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsor.control

## compare among treatments
dotchart(SES$SES_bsor, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_bsor ~ trmt))
with(SES, ad.test(SES_bsor))

####FS_bsor.lm <-lmer(SES_bsor ~ trmt+(trmt|site),data = SES)


FS_bsor.mod <- glm(SES_bsor ~ trmt , family = gaussian, data = SES)
summary(FS_bsor.mod)
qqnorm(resid(FS_bsor.mod))
qqline(resid(FS_bsor.mod))
plot(simulateResiduals(FS_bsor.mod))
densityPlot(rstudent(FS_bsor.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_bsor.mod)
influenceIndexPlot(FS_bsor.mod, vars = c("Cook"), id = list(n = 3))

summary(FS_bsor.mod)
Anova(FS_bsor.mod)
#No sig diff between treatments


### taxonomic diveristy - beta sim----
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

## compare among treatments
dotchart(SES$SES_bsim, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_bsim ~ trmt))
with(SES, ad.test(SES_bsim))

FS_bsim.mod <- glm(SES_bsim ~ trmt , family = gaussian, data = SES)
summary(FS_bsim.mod)
qqnorm(resid(FS_bsim.mod))
qqline(resid(FS_bsim.mod))
plot(simulateResiduals(FS_bsim.mod))
densityPlot(rstudent(FS_bsim.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_bsim.mod)
influenceIndexPlot(FS_bsim.mod, vars = c("Cook"), id = list(n = 3))

summary(FS_bsim.mod)
Anova(FS_bsim.mod)
#no sig difference


### taxonomic diveristy - beta sne----
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

## compare among treatments
dotchart(SES$SES_bsne, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_bsne ~ trmt))
with(SES, ad.test(SES_bsne))

FS_bsne.mod <- glm(SES_bsne ~ trmt , family = gaussian, data = SES)
summary(FS_bsne.mod)
qqnorm(resid(FS_bsne.mod))
qqline(resid(FS_bsne.mod))
plot(simulateResiduals(FS_bsne.mod))
densityPlot(rstudent(FS_bsne.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_bsne.mod)
influenceIndexPlot(FS_bsne.mod, vars = c("Cook"), id = list(n = 3))

summary(FS_bsne.mod)
Anova(FS_bsne.mod)
#no sig difference

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

FS_falpha.mod <- glm(SES_falpha ~ trmt , family = gaussian, data = SES)
summary(FS_falpha.mod)
qqnorm(resid(FS_falpha.mod))
qqline(resid(FS_falpha.mod))
plot(simulateResiduals(FS_falpha.mod))
densityPlot(rstudent(FS_falpha.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_falpha.mod)
influenceIndexPlot(FS_falpha.mod, vars = c("Cook"), id = list(n = 3))

summary(FS_falpha.mod)
Anova(FS_falpha.mod)

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

FS_fbsor.mod <- glm(SES_fbsor~ trmt , family = gaussian, data = SES)
summary(FS_fbsor.mod)
qqnorm(resid(FS_fbsor.mod))
qqline(resid(FS_fbsor.mod))
plot(simulateResiduals(FS_fbsor.mod))
densityPlot(rstudent(FS_fbsor.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_fbsor.mod)
influenceIndexPlot(FS_fbsor.mod, vars = c("Cook"), id = list(n = 3))

summary(FS_fbsor.mod)
Anova(FS_fbsor.mod)
#not significant by treatment

### functional beta diversity - beta sim----
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

## compare among treatments
dotchart(SES$SES_fbsim, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_fbsim ~ trmt))
with(SES, ad.test(SES_fbsim))

FS_fbsim.mod <- glm(SES_fbsim ~ trmt , family = gaussian, data = SES)
summary(FS_fbsim.mod)
qqnorm(resid(FS_fbsim.mod))
qqline(resid(FS_fbsim.mod))
plot(simulateResiduals(FS_fbsim.mod))
densityPlot(rstudent(FS_fbsim.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_fbsim.mod)
influenceIndexPlot(FS_fbsim.mod, vars = c("Cook"), id = list(n = 3))

summary(FS_fbsim.mod)
Anova(FS_fbsim.mod)


### functional beta diversity - beta sne----
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



## body length----
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

## compare among treatments
dotchart(SES$SES_bl, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_bl ~ trmt))
with(SES, ad.test(SES_bl))

FS_bl.mod <- glm(SES_bl ~ trmt , family = gaussian, data = SES)
summary(FS_bl.mod)
qqnorm(resid(FS_bl.mod))
qqline(resid(FS_bl.mod))
plot(simulateResiduals(FS_bl.mod))
densityPlot(rstudent(FS_bl.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_bl.mod)
influenceIndexPlot(FS_bl.mod, vars = c("Cook"), id = list(n = 3))

summary(FS_bl.mod)
Anova(FS_bl.mod)
#If significant run the below code and use the boxplot on figure code to make a graph. if not don't
emmeans(FS_bl.mod, pairwise ~ trmt)


## lecty---- 
###lec_0- Kleptoparasitic----
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

## compare among treatments
dotchart(SES$SES_lec_0, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_lec_0 ~ trmt))
with(SES, ad.test(SES_lec_0))

FS_lec_0.mod <- glm(SES_lec_0 ~ trmt , family = gaussian, data = SES)
summary(FS_lec_0.mod)
qqnorm(resid(FS_lec_0.mod))
qqline(resid(FS_lec_0.mod))
plot(simulateResiduals(FS_lec_0.mod))
densityPlot(rstudent(FS_lec_0.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_lec_0.mod)
influenceIndexPlot(FS_lec_0.mod, vars = c("Cook"), id = list(n = 3))

summary(FS_lec_0.mod)
Anova(FS_lec_0.mod)



### lec_1 - Generalist----
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

## compare among treatments
dotchart(SES$SES_lec_1, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_lec_1 ~ trmt))
with(SES, ad.test(SES_lec_1))

FS_lec_1.mod <- glm(SES_lec_1 ~ trmt , family = gaussian, data = SES)
summary(FS_lec_1.mod)
qqnorm(resid(FS_lec_1.mod))
qqline(resid(FS_lec_1.mod))
plot(simulateResiduals(FS_lec_1.mod))
densityPlot(rstudent(FS_lec_1.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_lec_1.mod)
influenceIndexPlot(FS_lec_1.mod, vars = c("Cook"), id = list(n = 3))

summary(FS_lec_1.mod)
Anova(FS_lec_1.mod)

### lec_2 - Specialist----
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


## compare among treatments
dotchart(SES$SES_lec_2, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_lec_2 ~ trmt))
with(SES, ad.test(SES_lec_2))

FS_lec_2.mod <- glm(SES_lec_2 ~ trmt , family = gaussian, data = SES)
summary(FS_lec_2.mod)
qqnorm(resid(FS_lec_2.mod))
qqline(resid(FS_lec_2.mod))
plot(simulateResiduals(FS_lec_2.mod))
densityPlot(rstudent(FS_lec_2.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_lec_2.mod)
influenceIndexPlot(FS_lec_2.mod, vars = c("Cook"), id = list(n = 3))
#There is an outlier detected- 47th st garden
FS_lec_2.mod.red <- update(FS_lec_2.mod, subset = -c(1))
summary(FS_lec_2.mod.red)
compareCoefs(FS_lec_2.mod, FS_lec_2.mod.red)
outlierTest(FS_lec_2.mod.red)
influenceIndexPlot(FS_lec_2.mod.red, vars = c("Cook"), id = list(n = 3))


summary(FS_lec_2.mod)
Anova(FS_lec_2.mod)
summary(FS_lec_2.mod.red)
Anova(FS_lec_2.mod.red)

#Not significant here 


##Nesting----
### nest_1 - Soil----
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

## compare among treatments
dotchart(SES$SES_nest_1, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_nest_1 ~ trmt))
with(SES, ad.test(SES_nest_1))

FS_nest_1.mod <- glm(SES_nest_1 ~ trmt , family = gaussian, data = SES)
summary(FS_nest_1.mod)
qqnorm(resid(FS_nest_1.mod))
qqline(resid(FS_nest_1.mod))
plot(simulateResiduals(FS_nest_1.mod))
densityPlot(rstudent(FS_nest_1.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_nest_1.mod)
influenceIndexPlot(FS_nest_1.mod, vars = c("Cook"), id = list(n = 3))
#There is an outlier detected- 47th st garden

summary(FS_nest_1.mod)
Anova(FS_nest_1.mod)
#The ANOVA said there was significance
emmeans(FS_nest_1.mod, pairwise ~ trmt)
#Significantly different- Farm looks like it has greater SES values than control


### nest_2 - Cavity----
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


## compare among treatments
dotchart(SES$SES_nest_2, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_nest_2 ~ trmt))
with(SES, ad.test(SES_nest_2))

FS_nest_2.mod <- glm(SES_nest_2 ~ trmt , family = gaussian, data = SES)
summary(FS_nest_2.mod)
qqnorm(resid(FS_nest_2.mod))
qqline(resid(FS_nest_2.mod))
plot(simulateResiduals(FS_nest_2.mod))
densityPlot(rstudent(FS_nest_2.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_nest_2.mod)
influenceIndexPlot(FS_nest_2.mod, vars = c("Cook"), id = list(n = 3))


summary(FS_nest_2.mod)
Anova(FS_nest_2.mod)
#The ANOVA said there was significance
emmeans(FS_nest_2.mod, pairwise ~ trmt)
#Significantly different- Farm looks like it has lower SES values than control
#THough both generally below 0


### nest_3 - Colony----
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

## compare among treatments
dotchart(SES$SES_nest_3, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_nest_3 ~ trmt))
with(SES, ad.test(SES_nest_3))

FS_nest_3.mod <- glm(SES_nest_3 ~ trmt , family = gaussian, data = SES)
summary(FS_nest_3.mod)
qqnorm(resid(FS_nest_3.mod))
qqline(resid(FS_nest_3.mod))
plot(simulateResiduals(FS_nest_3.mod))
densityPlot(rstudent(FS_nest_3.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_nest_3.mod)
influenceIndexPlot(FS_nest_3.mod, vars = c("Cook"), id = list(n = 3))


summary(FS_nest_3.mod)
Anova(FS_nest_3.mod)


### nest_4 - Pithy Stems----
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

## compare among treatments
dotchart(SES$SES_nest_4, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_nest_4 ~ trmt))
with(SES, ad.test(SES_nest_4))
#fails normality


FS_nest_4.mod <- glm(SES_nest_4 ~ trmt , family = gaussian, data = SES)
summary(FS_nest_4.mod)
qqnorm(resid(FS_nest_4.mod))
qqline(resid(FS_nest_4.mod))
plot(simulateResiduals(FS_nest_4.mod))
densityPlot(rstudent(FS_nest_4.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_nest_4.mod)
influenceIndexPlot(FS_nest_4.mod, vars = c("Cook"), id = list(n = 3))


summary(FS_nest_4.mod)
Anova(FS_nest_4.mod)
#no sig difference


### nest_5 - Wood----
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

## compare among treatments
dotchart(SES$SES_nest_5, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_nest_5 ~ trmt))
with(SES, ad.test(SES_nest_5))
#

FS_nest_5.mod <- glm(SES_nest_5 ~ trmt , family = gaussian, data = SES)
summary(FS_nest_5.mod)
qqnorm(resid(FS_nest_5.mod))
qqline(resid(FS_nest_5.mod))
plot(simulateResiduals(FS_nest_5.mod))
densityPlot(rstudent(FS_nest_5.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_nest_5.mod)
influenceIndexPlot(FS_nest_5.mod, vars = c("Cook"), id = list(n = 3))


summary(FS_nest_5.mod)
Anova(FS_nest_5.mod)
#no sig difference


##Sociality----
### soc_1 - Subsocial----
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


## compare among treatments
dotchart(SES$SES_soc_1, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_soc_1 ~ trmt))
with(SES, ad.test(SES_soc_1))
#

FS_soc_1.mod <- glm(SES_soc_1 ~ trmt , family = gaussian, data = SES)
summary(FS_soc_1.mod)
qqnorm(resid(FS_soc_1.mod))
qqline(resid(FS_soc_1.mod))
plot(simulateResiduals(FS_soc_1.mod))
densityPlot(rstudent(FS_soc_1.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_soc_1.mod)
influenceIndexPlot(FS_soc_1.mod, vars = c("Cook"), id = list(n = 3))


summary(FS_soc_1.mod)
Anova(FS_soc_1.mod)
#no sig difference



### soc_2 - Solitary----
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

## compare among treatments
dotchart(SES$SES_soc_2, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_soc_2 ~ trmt))
with(SES, ad.test(SES_soc_2))
#

FS_soc_2.mod <- glm(SES_soc_2 ~ trmt , family = gaussian, data = SES)
summary(FS_soc_2.mod)
qqnorm(resid(FS_soc_2.mod))
qqline(resid(FS_soc_2.mod))
plot(simulateResiduals(FS_soc_2.mod))
densityPlot(rstudent(FS_soc_2.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_soc_2.mod)
influenceIndexPlot(FS_soc_2.mod, vars = c("Cook"), id = list(n = 3))


summary(FS_soc_2.mod)
Anova(FS_soc_2.mod)
#no sig difference



### soc_3 - Eusocial----
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

## compare among treatments
dotchart(SES$SES_soc_3, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_soc_3 ~ trmt))
with(SES, ad.test(SES_soc_3))
#

FS_soc_3.mod <- glm(SES_soc_3 ~ trmt , family = gaussian, data = SES)
summary(FS_soc_3.mod)
qqnorm(resid(FS_soc_3.mod))
qqline(resid(FS_soc_3.mod))
plot(simulateResiduals(FS_soc_3.mod))
densityPlot(rstudent(FS_soc_3.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_soc_3.mod)
influenceIndexPlot(FS_soc_3.mod, vars = c("Cook"), id = list(n = 3))


summary(FS_soc_3.mod)
Anova(FS_soc_3.mod)
#no sig difference



### soc_4 - Parasitic----
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

## compare among treatments
dotchart(SES$SES_soc_4, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_soc_4 ~ trmt))
with(SES, ad.test(SES_soc_4))
#

FS_soc_4.mod <- glm(SES_soc_4 ~ trmt , family = gaussian, data = SES)
summary(FS_soc_4.mod)
qqnorm(resid(FS_soc_4.mod))
qqline(resid(FS_soc_4.mod))
plot(simulateResiduals(FS_soc_4.mod))
densityPlot(rstudent(FS_soc_4.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_soc_4.mod)
influenceIndexPlot(FS_soc_4.mod, vars = c("Cook"), id = list(n = 3))


summary(FS_soc_4.mod)
Anova(FS_soc_4.mod)
#no sig difference


##Origin----
### ori_0 - Native----
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

## compare among treatments
dotchart(SES$SES_ori_0, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_ori_0 ~ trmt))
with(SES, ad.test(SES_ori_0))
#

FS_ori_0.mod <- glm(SES_ori_0 ~ trmt , family = gaussian, data = SES)
summary(FS_ori_0.mod)
qqnorm(resid(FS_ori_0.mod))
qqline(resid(FS_ori_0.mod))
plot(simulateResiduals(FS_ori_0.mod))
densityPlot(rstudent(FS_ori_0.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_ori_0.mod)
influenceIndexPlot(FS_ori_0.mod, vars = c("Cook"), id = list(n = 3))


summary(FS_ori_0.mod)
Anova(FS_ori_0.mod)
#sig difference
emmeans(FS_ori_0.mod, pairwise ~ trmt)


### ori_1 - Exotic----
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


## compare among treatments
dotchart(SES$SES_ori_1, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_ori_1 ~ trmt))
with(SES, ad.test(SES_ori_1))
#

FS_ori_1.mod <- glm(SES_ori_1 ~ trmt , family = gaussian, data = SES)
summary(FS_ori_1.mod)
qqnorm(resid(FS_ori_1.mod))
qqline(resid(FS_ori_1.mod))
plot(simulateResiduals(FS_ori_1.mod))
densityPlot(rstudent(FS_ori_1.mod)) # check density estimate of the distribution of residuals
outlierTest(FS_ori_1.mod)
influenceIndexPlot(FS_ori_1.mod, vars = c("Cook"), id = list(n = 3))


summary(FS_ori_1.mod)
Anova(FS_ori_1.mod)
#sig difference
emmeans(FS_ori_1.mod, pairwise ~ trmt)



#Figures----

##Loading programs needed
if (!suppressWarnings(require(viridis))) install.packages("viridis")
citation("viridis")

if (!suppressWarnings(require(reshape2))) install.packages("reshape2")
citation("reshape2")


#step one- create a dataframe for all of the tax and funct diversity metrics I want to graph
SES_ALLdivcontrol <- as.data.frame(rbind( control$SES_fbsim, control$SES_fbsne, control$SES_fbsor, control$SES_bsim, control$SES_bsne, control$SES_bsor, control$SES_falpha))
str(SES_ALLdivcontrol)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdivcontrol <- data.frame(t(SES_ALLdivcontrol))
str(SES_ALLdivcontrol)
colnames(SES_ALLdivcontrol) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity","Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity", "Funct. Alpha" )


SES_ALLdivcontrol <- melt(SES_ALLdivcontrol)
colnames(SES_ALLdivcontrol) <- c("diversity","ses")



SES_ALLdivfarm <- as.data.frame(rbind( farm$SES_fbsim, farm$SES_fbsne, farm$SES_fbsor, farm$SES_bsim, farm$SES_bsne, farm$SES_bsor, farm$SES_falpha))
str(SES_ALLdivfarm)
SES_ALLdivfarm <- data.frame(t(SES_ALLdivfarm))
str(SES_ALLdivfarm)
colnames(SES_ALLdivfarm) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity","Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity", "Funct. Alpha")

SES_ALLdivfarm <- melt(SES_ALLdivfarm)
colnames(SES_ALLdivfarm) <- c("diversity","ses")


#Now I should be able to plot it in a nice graph
png("Figures/Figure4.png", width = 2500, height = 1200, pointsize = 20)
par(mfrow=c(1,2)) # indicates one row, two columns
par(mar = c(5,11,4,2)) # sets the margins around the figure

boxplot(ses ~ diversity, data = SES_ALLdivcontrol, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Vacant Lot Diversity", ylim=c(-21,18.5))
stripchart(ses ~ diversity, data = SES_ALLdivcontrol, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(15, 7.4, "A", pos = 4, font = 2, cex = 2.6)


boxplot(ses ~ diversity, data = SES_ALLdivfarm, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Farm Diversity", ylim=c(-21,18.5))
stripchart(ses ~ diversity, data = SES_ALLdivfarm, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(15, 7.4, "B", pos = 4, font = 2, cex = 2.6)

dev.off()


library(ggplot2)
library(cowplot)


#Figure 5----

#Length and Origin
SES_lengthandorigin.FC<- data.frame(
  trmt=rep(SES$trmt,3),
  variable=rep(c("Alien", "Native", "Body Length"), each = nrow(SES)),
  value = c(SES$SES_ori_1, SES$SES_ori_0, SES$SES_bl)
)
SES_lengthandorigin.FC$variable<-factor(SES_lengthandorigin.FC$variable, c("Alien", "Native", "Body Length"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots



fig5a<- ggplot(SES_lengthandorigin.FC, aes(x = variable, y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             shape = 19, size = 3, show.legend = TRUE) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=15)) +
  theme(axis.text.y = element_text(size=25))+
  ylim(c(-5, 5)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip() +
  theme(legend.text = element_text(size=18), legend.key.size = unit(2, 'cm'), legend.title = element_blank())




#Nesting traits

SES_nesting.FC<- data.frame(
  trmt=rep(SES$trmt,5),
  variable=rep(c("Wood","Pithy Stems", "Colony", "Cavity", "Soil"), each = nrow(SES)),
  value = c(SES$SES_nest_5, SES$SES_nest_4, SES$SES_nest_3, SES$SES_nest_2, SES$SES_nest_1)
)
SES_nesting.FC$variable<-factor(SES_nesting.FC$variable, c("Wood", "Pithy Stems","Colony", "Cavity", "Soil"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots


fig5b<- ggplot(SES_nesting.FC, aes(x = variable, y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             shape = 19, size = 3, show.legend = TRUE) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size=15)) +
  theme(axis.text.y = element_text(size=25))+
  ylim(c(-5, 5)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip()+
  theme(legend.text = element_text(size=18), legend.key.size = unit(2, 'cm'), legend.title = element_blank())



#Lecty
SES_lecty.FC<- data.frame(
  trmt=rep(SES$trmt,3),
  variable=rep(c("Specialist","Generalist", "Kleptoparasitic"), each = nrow(SES)),
  value = c(SES$SES_lec_2, SES$SES_lec_1, SES$SES_lec_0)
)
SES_lecty.FC$variable<-factor(SES_lecty.FC$variable, c("Specialist", "Generalist","Kleptoparasitic"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots


fig5c<- ggplot(SES_lecty.FC, aes(x = variable, y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             shape = 19, size = 3, show.legend = TRUE) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=15)) +
  theme(axis.text.y = element_text(size=25))+
  ylim(c(-5, 5)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip()+
  theme(legend.text = element_text(size=18), legend.key.size = unit(2, 'cm'), legend.title = element_blank())



#Sociality
SES_sociality.FC<- data.frame(
  trmt=rep(SES$trmt,4),
  variable=rep(c("Parasitic", "Eusocial","Subsocial", "Solitary"), each = nrow(SES)),
  value = c(SES$SES_soc_4, SES$SES_soc_3, SES$SES_soc_2, SES$SES_soc_1)
) #Note that we have to put our funct traits in the reverse order of how we want them to appear on our graph
SES_sociality.FC$variable<-factor(SES_sociality.FC$variable, c("Parasitic", "Eusocial","Subsocial", "Solitary"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots

fig5d<- ggplot(SES_sociality.FC, aes(x = variable, y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = .75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = .75), 
             shape = 19, size = 3) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") + 
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=15)) +
  theme(axis.text.y = element_text(size=25))+
  ylim(c(-5, 5)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip() +
  theme(legend.text = element_text(size=18), legend.key.size = unit(2, 'cm'), legend.title = element_blank())



fig5 <- plot_grid(fig5a,fig5b,fig5c,fig5d, labels= c('A','B','C','D'), label_size= 30)
png("Figures/Figure 5 all FS local sp traits.png", width = 1500, height = 1000, pointsize = 20)

fig5

dev.off()


treatmentcleanup.fs <- c(Farm="Urban Farms", Control="Vacant Lots")
SES$treatmentspelledout <- as.character(treatmentcleanup.fs[SES$trmt])

png("Figures/Figure 6.png", width = 1500, height = 1000, pointsize = 20)

par(mfrow=c(2,2)) # indicates three rows, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure

# Soil nesting
boxplot(SES_nest_1 ~ treatmentspelledout, data = SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-2,2.5), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Soil Nesting")
stripchart(SES_nest_1 ~ trmt, data = SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(2.38, 2.41, "A", pos = 4, font = 2, cex = 2)




# Cavity nesting
boxplot(SES_nest_2 ~ treatmentspelledout, data = SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-2,2.5), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Cavity Nesting")
stripchart(SES_nest_2 ~ trmt, data = SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(2.38, 2.41, "B", pos = 4, font = 2, cex = 2)

# native
boxplot(SES_ori_0 ~ treatmentspelledout, data = SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-2,2.5), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Native")
stripchart(SES_ori_0 ~ trmt, data = SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(2.38, 2.41, "C", pos = 4, font = 2, cex = 2)


# non-native
boxplot(SES_ori_1 ~ treatmentspelledout, data = SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        ylim = c(-2,2.5), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Non-Native")
stripchart(SES_ori_1 ~ trmt, data = SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(2.38, 2.41, "D", pos = 4, font = 2, cex = 2)

dev.off()
