###################################################################################
#
# Ohio bee data- MP
#
# Step 2: Urban Cleveland Pool to Local Greenspaces
# Comparing among Local habitats in MIchelle Pham's data
# 
#CWM and Functional Diversity - Null Models
# GLMMs
#
# 
#CA Shepard and KI Perry: 2 October 2023
###################################################################################

#Refer back to Final_local_analysis to get to the beginning of this document
#Everything before Line 83
t <- read.csv("traits_urbanpool.26.csv", row.names=1)
a <- read.csv("urbanpool.MP.26.csv", row.names=1)

aO<-a
a<-a[1:136]

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
write.csv(t, file = "Urban to Local_Nulls_MP/test.csv")
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
#checking that all sites have more than 4 species


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
wt <- c(0.38, 0.18, 0.11, 0.12, 0.21)

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
write.csv(nbl, file = "Urban to Local_Nulls_MP/nbl_MP.csv")
write.csv(nlec_0, file = "Urban to Local_Nulls_MP/nlec_0_MP.csv")
write.csv(nlec_1, file = "Urban to Local_Nulls_MP/nlec_1_MP.csv")
write.csv(nlec_2, file = "Urban to Local_Nulls_MP/nlec_2_MP.csv")
write.csv(nnest_1, file = "Urban to Local_Nulls_MP/nnest_1_MP.csv")
write.csv(nnest_2, file = "Urban to Local_Nulls_MP/nnest_2_MP.csv")
write.csv(nnest_3, file = "Urban to Local_Nulls_MP/nnest_3_MP.csv")
write.csv(nnest_4, file = "Urban to Local_Nulls_MP/nnest_4_MP.csv")
write.csv(nnest_5, file = "Urban to Local_Nulls_MP/nnest_5_MP.csv")
write.csv(nsoc_1, file = "Urban to Local_Nulls_MP/nsoc_1_MP.csv")
write.csv(nsoc_2, file = "Urban to Local_Nulls_MP/nsoc_2_MP.csv")
write.csv(nsoc_3, file = "Urban to Local_Nulls_MP/nsoc_3_MP.csv")
write.csv(nsoc_4, file = "Urban to Local_Nulls_MP/nsoc_4_MP.csv")
write.csv(nori_0, file = "Urban to Local_Nulls_MP/nori_0_MP.csv")
write.csv(nori_1, file = "Urban to Local_Nulls_MP/nori_1_MP.csv")

write.csv(nbsim, file = "Urban to Local_Nulls_MP/tbeta_sim_MP.csv")
write.csv(nbsne, file = "Urban to Local_Nulls_MP/tbeta_sne_MP.csv")
write.csv(nbsor, file = "Urban to Local_Nulls_MP/tbeta_sor_MP.csv")

write.csv(nfalpha, file = "Urban to Local_Nulls_MP/falpha_MP.csv")

write.csv(nfsim, file = "Urban to Local_Nulls_MP/fbeta_sim_MP.csv")
write.csv(nfsne, file = "Urban to Local_Nulls_MP/fbeta_sne_MP.csv")
write.csv(nfsor, file = "Urban to Local_Nulls_MP/fbeta_sor_MP.csv")

# load the output matrices

nbl <- read.csv("Urban to Local_Nulls_MP/nbl_MP.csv", row.names=1)
nlec_0 <- read.csv("Urban to Local_Nulls_MP/nlec_0_MP.csv", row.names=1)
nlec_1 <- read.csv("Urban to Local_Nulls_MP/nlec_1_MP.csv", row.names=1)
nlec_2 <- read.csv("Urban to Local_Nulls_MP/nlec_2_MP.csv", row.names=1)
nnest_1 <- read.csv("Urban to Local_Nulls_MP/nnest_1_MP.csv", row.names=1)
nnest_2 <- read.csv("Urban to Local_Nulls_MP/nnest_2_MP.csv", row.names=1)
nnest_3 <- read.csv("Urban to Local_Nulls_MP/nnest_3_MP.csv", row.names=1)
nnest_4 <- read.csv("Urban to Local_Nulls_MP/nnest_4_MP.csv", row.names=1)
nnest_5 <- read.csv("Urban to Local_Nulls_MP/nnest_5_MP.csv", row.names=1)
nsoc_1 <- read.csv("Urban to Local_Nulls_MP/nsoc_1_MP.csv", row.names=1)
nsoc_2 <- read.csv("Urban to Local_Nulls_MP/nsoc_2_MP.csv", row.names=1)
nsoc_3 <- read.csv("Urban to Local_Nulls_MP/nsoc_3_MP.csv", row.names=1)
nsoc_4 <- read.csv("Urban to Local_Nulls_MP/nsoc_4_MP.csv", row.names=1)
nori_0 <- read.csv("Urban to Local_Nulls_MP/nori_0_MP.csv", row.names=1)
nori_1 <- read.csv("Urban to Local_Nulls_MP/nori_1_MP.csv", row.names=1)

nfalpha <- read.csv("Urban to Local_Nulls_MP/falpha_MP.csv", row.names=1)

nbsim <- read.csv("Urban to Local_Nulls_MP/tbeta_sim_MP.csv", row.names=1)
nbsne <- read.csv("Urban to Local_Nulls_MP/tbeta_sne_MP.csv", row.names=1)
nbsor <- read.csv("Urban to Local_Nulls_MP/tbeta_sim_MP.csv", row.names=1)

nfsim <- read.csv("Urban to Local_Nulls_MP/fbeta_sim_MP.csv", row.names=1)
nfsne <- read.csv("Urban to Local_Nulls_MP/fbeta_sne_MP.csv", row.names=1)
nfsor <- read.csv("Urban to Local_Nulls_MP/fbeta_sor_MP.csv", row.names=1)

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

write.csv(SES.all, file = "Urban to Local_Nulls_MP/SES_Local_MP.csv")
#import the SES data
SES <- read.csv("Urban to Local_Nulls_MP/SES_Local_MP.csv", row.names = 1)

#Including treatment
a1<-aO
a1<-a1[-12,]
a1<-a1[-9,]
a1<-a1[-1,]

SES$trmt <- a1$trmt
str(SES)
SES

#Adding the random neighborhood variable
SES$neighd <- c("C", "DS", "F", "G", "H","SV","TM","SV","TM", "H", "SV", "BE", "C","DS","F","G")


## pull out data for each treatment

VL2019<- SES[which(SES$trmt== "VL"),]
str(VL2019)
UP2019<- SES[which(SES$trmt== "UP"),]

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
boxplot(SES$SES_bsor ~ SES$trmt, ylim = c(-25, 1))
abline(h = 0.0, col = "black", lwd = 3, lty=2)

# compare to null expectations by treatment
bsor.vl <- wilcox.test(VL2019$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsor.vl
bsor.up <- wilcox.test(UP2019$SES_bsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsor.up

## compare among treatments
dotchart(SES$SES_bsor, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_bsor ~ trmt))
with(SES, ad.test(SES_bsor))
#passed homogeneity test
#kruskal.test(SES_bsor ~ trmt, data = SES) #it was not significant here anyway

#Doing this to try and resolve an error I had with running the code below: install.packages("lme4", type = "source")
#library(lme4) 

mp_bsor.lm <- lmer(SES_bsor~trmt+(1|neighd), data = SES)
summary(mp_bsor.lm)
qqnorm(resid(mp_bsor.lm))
qqline(resid(mp_bsor.lm))
plot(simulateResiduals(mp_bsor.lm))
densityPlot(rstudent(mp_bsor.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_bsor.lm)
influenceIndexPlot(mp_bsor.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_bsor.lm)
Anova(mp_bsor.lm)
#not sig

#mp_bsor.mod <- glm(SES_bsor ~ trmt , family = gaussian, data = SES)
#summary(mp_bsor.mod)
#qqnorm(resid(mp_bsor.mod))
#qqline(resid(mp_bsor.mod))
#plot(simulateResiduals(mp_bsor.mod))
#densityPlot(rstudent(mp_bsor.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_bsor.mod)
#influenceIndexPlot(mp_bsor.mod, vars = c("Cook"), id = list(n = 3))

#summary(mp_bsor.mod)
#Anova(mp_bsor.mod)
#No sig diff between treatments



### taxonomic diveristy - beta sim----
hist(SES$SES_bsim)
plot(SES$SES_bsim)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
bsim.vl <- wilcox.test(VL2019$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsim.vl
bsim.up <- wilcox.test(UP2019$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsim.up

## compare among treatments
dotchart(SES$SES_bsim, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_bsim ~ trmt))
with(SES, ad.test(SES_bsim))
#passed homogeneity tests
#kruskal.test(SES_bsim ~ trmt, data = SES)#not significant anyway

mp_bsim.lm <- lmer(SES_bsim~trmt+(1|neighd), data = SES)
summary(mp_bsim.lm)
qqnorm(resid(mp_bsim.lm))
qqline(resid(mp_bsim.lm))
plot(simulateResiduals(mp_bsim.lm))
densityPlot(rstudent(mp_bsim.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_bsim.lm)
influenceIndexPlot(mp_bsim.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_bsim.lm)
Anova(mp_bsim.lm)
#not sig


#mp_bsim.mod <- glm(SES_bsim ~ trmt , family = gaussian, data = SES)
#summary(mp_bsim.mod)
#qqnorm(resid(mp_bsim.mod))
#qqline(resid(mp_bsim.mod))
#plot(simulateResiduals(mp_bsim.mod))
#densityPlot(rstudent(mp_bsim.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_bsim.mod)
#influenceIndexPlot(mp_bsim.mod, vars = c("Cook"), id = list(n = 3))

#summary(mp_bsim.mod)
#Anova(mp_bsim.mod)
#no sig difference


### taxonomic diveristy - beta sne----
hist(SES$SES_bsne)
plot(SES$SES_bsne)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bsne ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
bsne.vl <- wilcox.test(VL2019$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsne.vl
bsne.up <- wilcox.test(UP2019$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsne.up

## compare among treatments
dotchart(SES$SES_bsne, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_bsne ~ trmt))
with(SES, ad.test(SES_bsne))
#passed homogeneity test
#kruskal.test(SES_bsne ~ trmt, data = SES)
#not significant

mp_bsne.lm <- lmer(SES_bsne~trmt+(1|neighd), data = SES)
?isSingular
#okay based on that, I would say, not significant because an error came up indicating super low effects
#therefore: below is not needed, hastaging it

#summary(mp_bsne.lm)
#qqnorm(resid(mp_bsne.lm))
#qqline(resid(mp_bsne.lm))
#plot(simulateResiduals(mp_bsne.lm))
#densityPlot(rstudent(mp_bsne.lm)) # check density estimate of the distribution of residuals
#outlierTest(mp_bsne.lm)
#influenceIndexPlot(mp_bsne.lm, vars = c("Cook"), id = list(n = 3))

#summary(mp_bsne.lm)
Anova(mp_bsne.lm)
#not sig



#mp_bsne.mod <- glm(SES_bsne ~ trmt , family = gaussian, data = SES)
#summary(mp_bsne.mod)
#qqnorm(resid(mp_bsne.mod))
#qqline(resid(mp_bsne.mod))
#plot(simulateResiduals(mp_bsne.mod))
#densityPlot(rstudent(mp_bsne.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_bsne.mod)
#influenceIndexPlot(mp_bsne.mod, vars = c("Cook"), id = list(n = 3))

#summary(mp_bsne.mod)
#Anova(mp_bsne.mod)
#no sig difference


## functional alpha diversity----
hist(SES$SES_falpha)
plot(SES$SES_falpha)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_falpha ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)


## compare to null expectations by treatment
falpha.vl <- wilcox.test(VL2019$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
falpha.vl
falpha.up <- wilcox.test(UP2019$SES_falpha, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
falpha.up

## compare among treatments
dotchart(SES$SES_falpha, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_falpha ~ trmt))
with(SES, ad.test(SES_falpha))
#passed homogeneity test
#kruskal.test(SES_falpha ~ trmt, data = SES) #it was not significant here anyway

mp_falpha.lm <- lmer(SES_falpha~trmt+(1|neighd), data = SES)

summary(mp_falpha.lm)
qqnorm(resid(mp_falpha.lm))
qqline(resid(mp_falpha.lm))
plot(simulateResiduals(mp_falpha.lm))
densityPlot(rstudent(mp_falpha.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_falpha.lm)
influenceIndexPlot(mp_falpha.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_falpha.lm)
Anova(mp_falpha.lm)
#not sig


#mp_falpha.mod <- glm(SES_falpha ~ trmt , family = gaussian, data = SES)
#summary(mp_falpha.mod)
#qqnorm(resid(mp_falpha.mod))
#qqline(resid(mp_falpha.mod))
#plot(simulateResiduals(mp_falpha.mod))
#densityPlot(rstudent(mp_falpha.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_falpha.mod)
#influenceIndexPlot(mp_falpha.mod, vars = c("Cook"), id = list(n = 3))

#summary(mp_falpha.mod)
#Anova(mp_falpha.mod)


##Functional beta div----
## functional beta diversity - beta sor
hist(SES$SES_fbsor)
plot(SES$SES_fbsor, ylim = c(-0.5, 6.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsor ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
fbsor.VL <- wilcox.test(VL2019$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsor.VL
fbsor.UP <- wilcox.test(UP2019$SES_fbsor, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsor.UP

## compare among treatments
dotchart(SES$SES_fbsor, group = SES$trmt, pch = 19)
with(SES, bartlett.test(SES_fbsor ~ trmt))
with(SES, ad.test(SES_fbsor))
#passed homogeneity test
#kruskal.test(SES_fbsor ~ trmt, data = SES) #it was not significant here anyway

mp_fbsor.lm <- lmer(SES_fbsor~trmt+(1|neighd), data = SES)

summary(mp_falpha.lm)
qqnorm(resid(mp_falpha.lm))
qqline(resid(mp_falpha.lm))
plot(simulateResiduals(mp_falpha.lm))
densityPlot(rstudent(mp_falpha.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_falpha.lm)
influenceIndexPlot(mp_falpha.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_falpha.lm)
Anova(mp_falpha.lm)
#not sig


#mp_fbsor.mod <- glm(SES_fbsor~ trmt , family = gaussian, data = SES)
#summary(mp_fbsor.mod)
#qqnorm(resid(mp_fbsor.mod))
#qqline(resid(mp_fbsor.mod))
#plot(simulateResiduals(mp_fbsor.mod))
#densityPlot(rstudent(mp_fbsor.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_fbsor.mod)
#influenceIndexPlot(mp_fbsor.mod, vars = c("Cook"), id = list(n = 3))

#summary(mp_fbsor.mod)
#Anova(mp_fbsor.mod)
#not significant by treatment


### functional beta diversity - beta sim----
hist(SES$SES_fbsim)
plot(SES$SES_fbsim)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsim ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
fbsim.VL <- wilcox.test(VL2019$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsim.VL
fbsim.UP <- wilcox.test(UP2019$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsim.UP

## compare among treatments
dotchart(SES$SES_fbsim, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_fbsim ~ trmt))
with(SES, ad.test(SES_fbsim))

#kruskal.test(SES_fbsim ~ trmt, data = SES) 
#it was not significant 

mp_fbsim.lm <- lmer(SES_fbsim~trmt+(1|neighd), data = SES)
#boundary(singular)fit error
#assumed not significant because of above error

#summary(mp_fbsim.lm)
#qqnorm(resid(mp_fbsim.lm))
#qqline(resid(mp_fbsim.lm))
#plot(simulateResiduals(mp_fbsim.lm))
#densityPlot(rstudent(mp_fbsim.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_fbsim.lm)
influenceIndexPlot(mp_fbsim.lm, vars = c("Cook"), id = list(n = 3))

#summary(mp_fbsim.lm)
Anova(mp_fbsim.lm)
#not sig




### functional beta diversity - beta sne----
hist(SES$SES_fbsne)
plot(SES$SES_fbsne, ylim = c(-0.5, 6))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_fbsne ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
fbsne.VL <- wilcox.test(VL2019$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsne.VL
fbsne.UP <- wilcox.test(UP2019$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsne.UP


## compare among treatments
dotchart(SES$SES_fbsne, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_fbsne ~ trmt))
with(SES, ad.test(SES_fbsne))

#kruskal.test(SES_fbsne ~ trmt, data = SES) #it was not significant here anyway

mp_fbsne.lm <- lmer(SES_fbsne~trmt+(1|neighd), data = SES)
#
summary(mp_fbsne.lm)
qqnorm(resid(mp_fbsne.lm))
qqline(resid(mp_fbsne.lm))
plot(simulateResiduals(mp_fbsne.lm))
densityPlot(rstudent(mp_fbsne.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_fbsne.lm)
influenceIndexPlot(mp_fbsne.lm, vars = c("Cook"), id = list(n = 3))
#okay. outlier test had nothing BUT T7-h is OVER 1.2!!! on the cook plot
mp_fbsne.lm.red <- update(mp_fbsne.lm, subset = -c(10))
summary(mp_fbsne.lm.red)
compareCoefs(mp_fbsne.lm.red, mp_fbsne.lm)

summary(mp_fbsne.lm.red)
Anova(mp_fbsne.lm.red)
Anova(mp_fbsne.lm)
#not sig either way

#mp_fbsne.mod <- glm(SES_fbsne ~ trmt , family = gaussian, data = SES)
#summary(mp_fbsne.mod)
#qqnorm(resid(mp_fbsne.mod))
#qqline(resid(mp_fbsne.mod))
#plot(simulateResiduals(mp_fbsne.mod))
#densityPlot(rstudent(mp_fbsne.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_fbsne.mod)
#influenceIndexPlot(mp_fbsne.mod, vars = c("Cook"), id = list(n = 3))

#summary(mp_fbsne.mod)
#Anova(mp_fbsne.mod)
#not significant



## body length----
hist(SES$SES_bl)
plot(SES$SES_bl, pch = 19, cex = 1.5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_bl ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
bl.VL <- wilcox.test(VL2019$SES_bl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bl.VL
bl.UP <- wilcox.test(UP2019$SES_bl, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bl.UP

## compare among treatments
dotchart(SES$SES_bl, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_bl ~ trmt))
with(SES, ad.test(SES_bl))
#passed homogeneity test
#kruskal.test(SES_bl ~ trmt, data = SES) #it was not significant here anyway

mp_bl.lm <- lmer(SES_bl~trmt+(1|neighd), data = SES)
#
summary(mp_bl.lm)
qqnorm(resid(mp_bl.lm))
qqline(resid(mp_bl.lm))
plot(simulateResiduals(mp_bl.lm))
densityPlot(rstudent(mp_bl.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_bl.lm)
influenceIndexPlot(mp_bl.lm, vars = c("Cook"), id = list(n = 3))
#WOW T8-BE has a cooks distance of over 2.5??? but it ISNT an outlier
mp_bl.lm.red <- update(mp_bl.lm, subset = -c(12))
influenceIndexPlot(mp_bl.lm.red, vars = c("Cook"), id = list(n = 3))

summary(mp_bl.lm.red)
compareCoefs(mp_bl.lm.red, mp_bl.lm)



summary(mp_bl.lm)
Anova(mp_bl.lm.red)
Anova(mp_bl.lm)
#not sig either way


#mp_bl.mod <- glm(SES_bl ~ trmt , family = gaussian, data = SES)
##summary(mp_bl.mod)
#qqnorm(resid(mp_bl.mod))
#qqline(resid(mp_bl.mod))
#plot(simulateResiduals(mp_bl.mod))
#densityPlot(rstudent(mp_bl.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_bl.mod)
#influenceIndexPlot(mp_bl.mod, vars = c("Cook"), id = list(n = 3))

#summary(mp_bl.mod)
#Anova(mp_bl.mod)
#Not significant


## lecty---- 
###lec_0- Kleptoparasitic----
hist(SES$SES_lec_0)
plot(SES$SES_lec_0, ylim = c(-3.5, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_lec_0 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
lec_0.VL <- wilcox.test(VL2019$SES_lec_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_0.VL
lec_0.UP <- wilcox.test(UP2019$SES_lec_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_0.UP

## compare among treatments
dotchart(SES$SES_lec_0, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_lec_0 ~ trmt))
with(SES, ad.test(SES_lec_0))
#passed homogeneity test
#kruskal.test(SES_lec_0 ~ trmt, data = SES) 

mp_lec_0.lm <- lmer(SES_lec_0~trmt+(1|neighd), data = SES)

summary(mp_lec_0.lm)
qqnorm(resid(mp_lec_0.lm))
qqline(resid(mp_lec_0.lm))
plot(simulateResiduals(mp_lec_0.lm))
densityPlot(rstudent(mp_lec_0.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_lec_0.lm)
influenceIndexPlot(mp_lec_0.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_lec_0.lm)
Anova(mp_lec_0.lm)
#not sig


#mp_lec_0.mod <- glm(SES_lec_0 ~ trmt , family = gaussian, data = SES)
#summary(mp_lec_0.mod)
#qqnorm(resid(mp_lec_0.mod))
#qqline(resid(mp_lec_0.mod))
#plot(simulateResiduals(mp_lec_0.mod))
#densityPlot(rstudent(mp_lec_0.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_lec_0.mod)
#influenceIndexPlot(mp_lec_0.mod, vars = c("Cook"), id = list(n = 3))
#okay there are definitely 2 outliers that seem to have some influence
#mp_lec_0.mod.red <- update(mp_lec_0.mod, subset = -c(3,11))
#summary(mp_lec_0.mod.red)
#compareCoefs(mp_lec_0.mod.red, mp_lec_0.mod)
#changed quite a bit

#summary(mp_lec_0.mod)
#Anova(mp_lec_0.mod)
#summary(mp_lec_0.mod.red)
#Anova(mp_lec_0.mod.red)

#qqnorm(resid(mp_lec_0.mod.red))
#qqline(resid(mp_lec_0.mod.red))
#plot(simulateResiduals(mp_lec_0.mod.red))
#densityPlot(rstudent(mp_lec_0.mod.red)) # check density estimate of the distribution of residuals
#outlierTest(mp_lec_0.mod.red)
#influenceIndexPlot(mp_lec_0.mod.red, vars = c("Cook"), id = list(n = 3))

 


### lec_1 - Generalist----
hist(SES$SES_lec_1)
SES$SES_lec_1
plot(SES$SES_lec_1)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_lec_1 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
lec_1.VL <- wilcox.test(VL2019$SES_lec_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_1.VL
lec_1.UP <- wilcox.test(UP2019$SES_lec_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_1.UP

## compare among treatments
dotchart(SES$SES_lec_1, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_lec_1 ~ trmt))
with(SES, ad.test(SES_lec_1))
#passed homogeneity test
#kruskal.test(SES_lec_1 ~ trmt, data = SES) #it was not significant here anyway

mp_lec_1.lm <- lmer(SES_lec_1~trmt+(1|neighd), data = SES)
#
summary(mp_lec_1.lm)
qqnorm(resid(mp_lec_1.lm))
qqline(resid(mp_lec_1.lm))
plot(simulateResiduals(mp_lec_1.lm))
densityPlot(rstudent(mp_lec_1.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_lec_1.lm)
influenceIndexPlot(mp_lec_1.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_lec_1.lm)
Anova(mp_lec_1.lm)
#not sig



#mp_lec_1.mod <- glm(SES_lec_1 ~ trmt , family = gaussian, data = SES)
#summary(mp_lec_1.mod)
#qqnorm(resid(mp_lec_1.mod))
#qqline(resid(mp_lec_1.mod))
#plot(simulateResiduals(mp_lec_1.mod))
#densityPlot(rstudent(mp_lec_1.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_lec_1.mod)
#influenceIndexPlot(mp_lec_1.mod, vars = c("Cook"), id = list(n = 3))

#summary(mp_lec_1.mod)
#Anova(mp_lec_1.mod)
#not significant



### lec_2 - Specialist----
hist(SES$SES_lec_2)
SES$SES_lec_2
plot(SES$SES_lec_2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_lec_2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
lec_2.VL <- wilcox.test(VL2019$SES_lec_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_2.VL
lec_2.UP <- wilcox.test(UP2019$SES_lec_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_2.UP


## compare among treatments
dotchart(SES$SES_lec_2, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_lec_2 ~ trmt))
with(SES, ad.test(SES_lec_2))
#passed homogeneity test
#kruskal.test(SES_lec_2 ~ trmt, data = SES) #it was not significant here anyway

mp_lec_2.lm <- lmer(SES_lec_2~trmt+(1|neighd), data = SES)
#boundary(signular)fit
#summary(mp_lec_2.lm)
#qqnorm(resid(mp_lec_2.lm))
#qqline(resid(mp_lec_2.lm))
#plot(simulateResiduals(mp_lec_2.lm))
#densityPlot(rstudent(mp_lec_2.lm)) # check density estimate of the distribution of residuals
#outlierTest(mp_lec_2.lm)
#influenceIndexPlot(mp_lec_2.lm, vars = c("Cook"), id = list(n = 3))

#summary(mp_lec_2.lm)
Anova(mp_lec_2.lm)
#not sig

#mp_lec_2.mod <- glm(SES_lec_2 ~ trmt , family = gaussian, data = SES)
#summary(mp_lec_2.mod)
#qqnorm(resid(mp_lec_2.mod))
#qqline(resid(mp_lec_2.mod))
#plot(simulateResiduals(mp_lec_2.mod))
#densityPlot(rstudent(mp_lec_2.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_lec_2.mod)
#influenceIndexPlot(mp_lec_2.mod, vars = c("Cook"), id = list(n = 3))
#No sig outlier

#summary(mp_lec_2.mod)
#Anova(mp_lec_2.mod)
#not significant


##Nesting----
### nest_1 - Soil----
hist(SES$SES_nest_1)
plot(SES$SES_nest_1)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest_1 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
nest_1.VL <- wilcox.test(VL2019$SES_nest_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_1.VL
nest_1.UP <- wilcox.test(UP2019$SES_nest_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_1.UP

## compare among treatments
dotchart(SES$SES_nest_1, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_nest_1 ~ trmt))
with(SES, ad.test(SES_nest_1))
#passed homogeneity test
#kruskal.test(SES_nest_1 ~ trmt, data = SES) #it was not significant- 0.0502

mp_nest_1.lm <- lmer(SES_nest_1~trmt+(1|neighd), data = SES)

summary(mp_nest_1.lm)
qqnorm(resid(mp_nest_1.lm))
qqline(resid(mp_nest_1.lm))
plot(simulateResiduals(mp_nest_1.lm))
densityPlot(rstudent(mp_nest_1.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_nest_1.lm)
influenceIndexPlot(mp_nest_1.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_nest_1.lm)
Anova(mp_nest_1.lm)
#Okay that is NOt significant!


#mp_nest_1.mod <- glm(SES_nest_1 ~ trmt , family = gaussian, data = SES)
#summary(mp_nest_1.mod)
#qqnorm(resid(mp_nest_1.mod))
#qqline(resid(mp_nest_1.mod))
#plot(simulateResiduals(mp_nest_1.mod))
#densityPlot(rstudent(mp_nest_1.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_nest_1.mod)
#influenceIndexPlot(mp_nest_1.mod, vars = c("Cook"), id = list(n = 3))


#summary(mp_nest_1.mod)
#Anova(mp_nest_1.mod)
#significant difference with more soil nesting bees in vacant ltos than urban prairies


### nest_2 - Cavity----
hist(SES$SES_nest_2)
plot(SES$SES_nest_2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest_2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
nest_2.VL <- wilcox.test(VL2019$SES_nest_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_2.VL
nest_2.UP <- wilcox.test(UP2019$SES_nest_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_2.UP


## compare among treatments
dotchart(SES$SES_nest_2, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_nest_2 ~ trmt))
with(SES, ad.test(SES_nest_2))
#passed normality  test
kruskal.test(SES_nest_2 ~ trmt, data = SES) #it was not significant here anyway

mp_nest_2.lm <- lmer(SES_nest_2~trmt+(1|neighd), data = SES)
#
summary(mp_nest_2.lm)
qqnorm(resid(mp_nest_2.lm))
qqline(resid(mp_nest_2.lm))
plot(simulateResiduals(mp_nest_2.lm))
densityPlot(rstudent(mp_nest_2.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_nest_2.lm)
influenceIndexPlot(mp_nest_2.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_nest_2.lm)
Anova(mp_nest_2.lm)
#not significant


#mp_nest_2.mod <- glm(SES_nest_2 ~ trmt , family = gaussian, data = SES)
#summary(mp_nest_2.mod)
#qqnorm(resid(mp_nest_2.mod))
#qqline(resid(mp_nest_2.mod))
#plot(simulateResiduals(mp_nest_2.mod))
#densityPlot(rstudent(mp_nest_2.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_nest_2.mod)
#influenceIndexPlot(mp_nest_2.mod, vars = c("Cook"), id = list(n = 3))


#summary(mp_nest_2.mod)
#Anova(mp_nest_2.mod)
#no significance


### nest_3 - Colony----
hist(SES$SES_nest_3)
plot(SES$SES_nest_3)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest_3 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
nest_3.VL <- wilcox.test(VL2019$SES_nest_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_3.VL
nest_3.UP <- wilcox.test(UP2019$SES_nest_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_3.UP

## compare among treatments
dotchart(SES$SES_nest_3, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_nest_3 ~ trmt))
with(SES, ad.test(SES_nest_3))
#passed homogeneity test
#kruskal.test(SES_nest_3 ~ trmt, data = SES) #it was not significant here anyway

mp_nest_3.lm <- lmer(SES_nest_3~trmt+(1|neighd), data = SES)
#
summary(mp_nest_3.lm)
qqnorm(resid(mp_nest_3.lm))
qqline(resid(mp_nest_3.lm))
plot(simulateResiduals(mp_nest_3.lm))
densityPlot(rstudent(mp_nest_3.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_nest_3.lm)
influenceIndexPlot(mp_nest_3.lm, vars = c("Cook"), id = list(n = 3))
#OKAY T8-BE has an influence of over 4
mp_nest_3.lm.red <- update(mp_nest_3.lm, subset = -c(12))
influenceIndexPlot(mp_nest_3.lm.red, vars = c("Cook"), id = list(n = 3))

summary(mp_nest_3.lm.red)
compareCoefs(mp_nest_3.lm.red, mp_nest_3.lm)



summary(mp_nest_3.lm)
Anova(mp_nest_3.lm.red)
Anova(mp_nest_3.lm)

#not sig either way


#mp_nest_3.mod <- glm(SES_nest_3 ~ trmt , family = gaussian, data = SES)
#summary(mp_nest_3.mod)
#qqnorm(resid(mp_nest_3.mod))
#qqline(resid(mp_nest_3.mod))
#plot(simulateResiduals(mp_nest_3.mod))
#densityPlot(rstudent(mp_nest_3.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_nest_3.mod)
#influenceIndexPlot(mp_nest_3.mod, vars = c("Cook"), id = list(n = 3))


#summary(mp_nest_3.mod)
#Anova(mp_nest_3.mod)



### nest_4 - Pithy Stems----
hist(SES$SES_nest_4)
plot(SES$SES_nest_4)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest_4 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
nest_4.VL <- wilcox.test(VL2019$SES_nest_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_4.VL
nest_4.UP <- wilcox.test(UP2019$SES_nest_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_4.UP

## compare among treatments
dotchart(SES$SES_nest_4, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_nest_4 ~ trmt))
with(SES, ad.test(SES_nest_4))
# passed normality test
kruskal.test(SES_nest_4 ~ trmt, data = SES) #it was not significant here 

mp_nest_4.lm <- lmer(SES_nest_4~trmt+(1|neighd), data = SES)
#
summary(mp_nest_4.lm)
qqnorm(resid(mp_nest_4.lm))
qqline(resid(mp_nest_4.lm))
plot(simulateResiduals(mp_nest_4.lm))
densityPlot(rstudent(mp_nest_4.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_nest_4.lm)
influenceIndexPlot(mp_nest_4.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_nest_4.lm)
Anova(mp_nest_4.lm)
#not sig


#mp_nest_4.mod <- glm(SES_nest_4 ~ trmt , family = gaussian, data = SES)
#summary(mp_nest_4.mod)
#qqnorm(resid(mp_nest_4.mod))
#qqline(resid(mp_nest_4.mod))
#plot(simulateResiduals(mp_nest_4.mod))
#densityPlot(rstudent(mp_nest_4.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_nest_4.mod)
##There is an outlier- influence plot looks like  T8-D
#influenceIndexPlot(mp_nest_4.mod, vars = c("Cook"), id = list(n = 3))
#mp_nest_4.mod.red <- update(mp_nest_4.mod, subset = -c(14))
#summary(mp_nest_4.mod.red)
#compareCoefs(mp_nest_4.mod.red, mp_nest_4.mod)


#summary(mp_nest_4.mod)
#Anova(mp_nest_4.mod)
#Anova(mp_nest_4.mod.red)
#no sig difference


### nest_5 - Wood----
hist(SES$SES_nest_5)
plot(SES$SES_nest_5)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_nest_5 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
nest_5.VL <- wilcox.test(VL2019$SES_nest_5, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_5.VL
nest_5.UP <- wilcox.test(UP2019$SES_nest_5, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_5.UP

## compare among treatments
dotchart(SES$SES_nest_5, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_nest_5 ~ trmt))
with(SES, ad.test(SES_nest_5))
#passed homogeneity test
#kruskal.test(SES_nest_5 ~ trmt, data = SES) #it was not significant here anyway

mp_nest_5.lm <- lmer(SES_nest_5~trmt+(1|neighd), data = SES)
#
summary(mp_nest_5.lm)
qqnorm(resid(mp_nest_5.lm))
qqline(resid(mp_nest_5.lm))
plot(simulateResiduals(mp_nest_5.lm))
densityPlot(rstudent(mp_nest_5.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_nest_5.lm)
influenceIndexPlot(mp_nest_5.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_nest_5.lm)
Anova(mp_nest_5.lm)
#not significant


#mp_nest_5.mod <- glm(SES_nest_5 ~ trmt , family = gaussian, data = SES)
#summary(mp_nest_5.mod)
#qqnorm(resid(mp_nest_5.mod))
#qqline(resid(mp_nest_5.mod))
#plot(simulateResiduals(mp_nest_5.mod))
#densityPlot(rstudent(mp_nest_5.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_nest_5.mod)
#influenceIndexPlot(mp_nest_5.mod, vars = c("Cook"), id = list(n = 3))


#summary(mp_nest_5.mod)
#Anova(mp_nest_5.mod)
#no sig difference


##Sociality----
### soc_1 - Subsocial----
hist(SES$SES_soc_1)
plot(SES$SES_soc_1)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_soc_1 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
soc_1.VL <- wilcox.test(VL2019$SES_soc_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_1.VL
soc_1.UP <- wilcox.test(UP2019$SES_soc_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_1.UP


## compare among treatments
dotchart(SES$SES_soc_1, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_soc_1 ~ trmt))
with(SES, ad.test(SES_soc_1))
# did not pass homogeneity test
kruskal.test(SES_soc_1 ~ trmt, data = SES) #it was not significant here 

mp_soc_1.lm <- lmer(SES_soc_1~trmt+(1|neighd), data = SES)
#boundary(signular)fit
summary(mp_soc_1.lm)
qqnorm(resid(mp_soc_1.lm))
qqline(resid(mp_soc_1.lm))
plot(simulateResiduals(mp_soc_1.lm))
densityPlot(rstudent(mp_soc_1.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_soc_1.lm)
influenceIndexPlot(mp_soc_1.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_soc_1.lm)
Anova(mp_soc_1.lm)
#not significant


#mp_soc_1.mod <- glm(SES_soc_1 ~ trmt , family = gaussian, data = SES)
#summary(mp_soc_1.mod)
#qqnorm(resid(mp_soc_1.mod))
#qqline(resid(mp_soc_1.mod))
#plot(simulateResiduals(mp_soc_1.mod))
#densityPlot(rstudent(mp_soc_1.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_soc_1.mod)
#influenceIndexPlot(mp_soc_1.mod, vars = c("Cook"), id = list(n = 3))

#summary(mp_soc_1.mod)
#Anova(mp_soc_1.mod)
#no sig difference


### soc_2 - Solitary----
hist(SES$SES_soc_2)
plot(SES$SES_soc_2)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_soc_2 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
soc_2.VL <- wilcox.test(VL2019$SES_soc_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_2.VL
soc_2.UP <- wilcox.test(UP2019$SES_soc_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_2.UP

## compare among treatments
dotchart(SES$SES_soc_2, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_soc_2 ~ trmt))
with(SES, ad.test(SES_soc_2))
#failed homogeneity test
kruskal.test(SES_soc_2 ~ trmt, data = SES) #almost significant here 
oneway.test(SES_soc_2 ~ trmt, data = SES, var.equal = FALSE)#Not significant


mp_soc_2.lm <- lmer(SES_soc_2~trmt+(1|neighd), data = SES)
#boundary(signular)fit
summary(mp_soc_2.lm)
qqnorm(resid(mp_soc_2.lm))
qqline(resid(mp_soc_2.lm))
plot(simulateResiduals(mp_soc_2.lm))
densityPlot(rstudent(mp_soc_2.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_soc_2.lm)
influenceIndexPlot(mp_soc_2.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_soc_2.lm)
Anova(mp_soc_2.lm)
#not significant


#mp_soc_2.mod <- glm(SES_soc_2 ~ trmt , family = gaussian, data = SES)
#summary(mp_soc_2.mod)
#qqnorm(resid(mp_soc_2.mod))
#qqline(resid(mp_soc_2.mod))
#plot(simulateResiduals(mp_soc_2.mod))
#densityPlot(rstudent(mp_soc_2.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_soc_2.mod)
#influenceIndexPlot(mp_soc_2.mod, vars = c("Cook"), id = list(n = 3))


#summary(mp_soc_2.mod)
#Anova(mp_soc_2.mod)
#no sig difference



### soc_3 - Eusocial----
hist(SES$SES_soc_3)
plot(SES$SES_soc_3)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_soc_3 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
soc_3.VL <- wilcox.test(VL2019$SES_soc_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_3.VL
soc_3.UP <- wilcox.test(UP2019$SES_soc_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_3.UP

## compare among treatments
dotchart(SES$SES_soc_3, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_soc_3 ~ trmt))
with(SES, ad.test(SES_soc_3))
#passed homogeneity test
#kruskal.test(SES_soc_3 ~ trmt, data = SES) #it was not significant here anyway

mp_soc_3.lm <- lmer(SES_soc_3~trmt+(1|neighd), data = SES)
#
summary(mp_soc_3.lm)
qqnorm(resid(mp_soc_3.lm))
qqline(resid(mp_soc_3.lm))
plot(simulateResiduals(mp_soc_3.lm))
densityPlot(rstudent(mp_soc_3.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_soc_3.lm)
influenceIndexPlot(mp_soc_3.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_soc_3.lm)
Anova(mp_soc_3.lm)
#not significant


#mp_soc_3.mod <- glm(SES_soc_3 ~ trmt , family = gaussian, data = SES)
#summary(mp_soc_3.mod)
#qqnorm(resid(mp_soc_3.mod))
#qqline(resid(mp_soc_3.mod))
#plot(simulateResiduals(mp_soc_3.mod))
#densityPlot(rstudent(mp_soc_3.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_soc_3.mod)
#influenceIndexPlot(mp_soc_3.mod, vars = c("Cook"), id = list(n = 3))


#summary(mp_soc_3.mod)
#Anova(mp_soc_3.mod)
#no sig difference



### soc_4 - Parasitic----
hist(SES$SES_soc_4)
plot(SES$SES_soc_4, ylim = c(-3.5, 0.5))
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_soc_4 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
soc_4.VL <- wilcox.test(VL2019$SES_soc_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_4.VL
soc_4.UP <- wilcox.test(UP2019$SES_soc_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_4.UP

## compare among treatments
dotchart(SES$SES_soc_4, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_soc_4 ~ trmt))
with(SES, ad.test(SES_soc_4))
#passed homogeneity test
#kruskal.test(SES_soc_4 ~ trmt, data = SES) #it was significant  anyway

mp_soc_4.lm <- lmer(SES_soc_4~trmt+(1|neighd), data = SES)
#
summary(mp_soc_4.lm)
qqnorm(resid(mp_soc_4.lm))
qqline(resid(mp_soc_4.lm))
plot(simulateResiduals(mp_soc_4.lm))
densityPlot(rstudent(mp_soc_4.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_soc_4.lm)
influenceIndexPlot(mp_soc_4.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_soc_4.lm)
Anova(mp_soc_4.lm)
#NO LONGER a Significant difference! (UPDATE WITH 2025 RUN)


#mp_soc_4.mod <- glm(SES_soc_4 ~ trmt , family = gaussian, data = SES)
#summary(mp_soc_4.mod)
#qqnorm(resid(mp_soc_4.mod))
#qqline(resid(mp_soc_4.mod))
#plot(simulateResiduals(mp_soc_4.mod))
#densityPlot(rstudent(mp_soc_4.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_soc_4.mod)
#influenceIndexPlot(mp_soc_4.mod, vars = c("Cook"), id = list(n = 3))
#might be two outliers, but one from each treatment

#summary(mp_soc_4.mod)
#Anova(mp_soc_4.mod)
# sig difference



##Origin----
### ori_0 - Native----
hist(SES$SES_ori_0)
plot(SES$SES_ori_0)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_ori_0 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
ori_0.VL <- wilcox.test(VL2019$SES_ori_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_0.VL
ori_0.UP <- wilcox.test(UP2019$SES_ori_0, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_0.UP

## compare among treatments
dotchart(SES$SES_ori_0, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_ori_0 ~ trmt))
with(SES, ad.test(SES_ori_0))
#passed homogeneity test
#kruskal.test(SES_ori_0 ~ trmt, data = SES) #it was not significant here anyway

mp_ori_0.lm <- lmer(SES_ori_0~trmt+(1|neighd), data = SES)
#boundary(signular)fit
summary(mp_ori_0.lm)
qqnorm(resid(mp_ori_0.lm))
qqline(resid(mp_ori_0.lm))
plot(simulateResiduals(mp_ori_0.lm))
densityPlot(rstudent(mp_ori_0.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_ori_0.lm)
influenceIndexPlot(mp_ori_0.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_ori_0.lm)
Anova(mp_ori_0.lm)
#not significant


#mp_ori_0.mod <- glm(SES_ori_0 ~ trmt , family = gaussian, data = SES)
#summary(mp_ori_0.mod)
#qqnorm(resid(mp_ori_0.mod))
#qqline(resid(mp_ori_0.mod))
#plot(simulateResiduals(mp_ori_0.mod))
#densityPlot(rstudent(mp_ori_0.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_ori_0.mod)
#influenceIndexPlot(mp_ori_0.mod, vars = c("Cook"), id = list(n = 3))


#summary(mp_ori_0.mod)
#Anova(mp_ori_0.mod)
#no sig difference



### ori_2 - Exotic----
hist(SES$SES_ori_1)
plot(SES$SES_ori_1)
abline(h = 0.0, col = "black", lwd = 3, lty=2)
boxplot(SES$SES_ori_1 ~ SES$trmt)
abline(h = 0.0, col = "black", lwd = 3, lty=2)

## compare to null expectations by treatment
ori_1.VL <- wilcox.test(VL2019$SES_ori_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_1.VL
ori_1.UP <- wilcox.test(UP2019$SES_ori_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
ori_1.UP


## compare among treatments
dotchart(SES$SES_ori_1, group = SES$trmt, pch = 19)

with(SES, bartlett.test(SES_ori_1 ~ trmt))
with(SES, ad.test(SES_ori_1))
#passed homogeneity test
#kruskal.test(SES_ori_1 ~ trmt, data = SES) #it was not significant here anyway

mp_ori_1.lm <- lmer(SES_ori_1~trmt+(1|neighd), data = SES)
#boundary(signular)fit
summary(mp_ori_1.lm)
qqnorm(resid(mp_ori_1.lm))
qqline(resid(mp_ori_1.lm))
plot(simulateResiduals(mp_ori_1.lm))
densityPlot(rstudent(mp_ori_1.lm)) # check density estimate of the distribution of residuals
outlierTest(mp_ori_1.lm)
influenceIndexPlot(mp_ori_1.lm, vars = c("Cook"), id = list(n = 3))

summary(mp_ori_1.lm)
Anova(mp_ori_1.lm)
#not significant


#mp_ori_1.mod <- glm(SES_ori_1 ~ trmt , family = gaussian, data = SES)
#summary(mp_ori_1.mod)
#qqnorm(resid(mp_ori_1.mod))
##qqline(resid(mp_ori_1.mod))
#plot(simulateResiduals(mp_ori_1.mod))
#densityPlot(rstudent(mp_ori_1.mod)) # check density estimate of the distribution of residuals
#outlierTest(mp_ori_1.mod)
#influenceIndexPlot(mp_ori_1.mod, vars = c("Cook"), id = list(n = 3))


#summary(mp_ori_1.mod)
#Anova(mp_ori_1.mod)
#no sig dif




#Figures----

##Loading programs needed
if (!suppressWarnings(require(viridis))) install.packages("viridis")
citation("viridis")

if (!suppressWarnings(require(reshape2))) install.packages("reshape2")
citation("reshape2")

#step one- create a dataframe for all of the tax and funct diversity metrics I want to graph
SES_ALLdivvl <- as.data.frame(rbind( VL2019$SES_fbsim, VL2019$SES_fbsne, VL2019$SES_fbsor, VL2019$SES_bsim, VL2019$SES_bsne, VL2019$SES_bsor, VL2019$SES_falpha))
str(SES_ALLdivvl)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdivvl <- data.frame(t(SES_ALLdivvl))
str(SES_ALLdivvl)
colnames(SES_ALLdivvl) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity","Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity", "Funct. Alpha" )


SES_ALLdivvl <- melt(SES_ALLdivvl)
colnames(SES_ALLdivvl) <- c("diversity","ses")



SES_ALLdivup <- as.data.frame(rbind( UP2019$SES_fbsim, UP2019$SES_fbsne, UP2019$SES_fbsor, UP2019$SES_bsim, UP2019$SES_bsne, UP2019$SES_bsor, UP2019$SES_falpha))
str(SES_ALLdivup)
SES_ALLdivup <- data.frame(t(SES_ALLdivup))
str(SES_ALLdivup)
colnames(SES_ALLdivup) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity","Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity", "Funct. Alpha")

SES_ALLdivup <- melt(SES_ALLdivup)
colnames(SES_ALLdivup) <- c("diversity","ses")


#Now I should be able to plot it in a nice graph
png("Figures/Figure 10.png", width = 2500, height = 1200, pointsize = 20)
par(mfrow=c(1,2)) # indicates one row, two columns
par(mar = c(5,11,4,2)) # sets the margins around the figure

boxplot(ses ~ diversity, data = SES_ALLdivvl, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "2019 Vacant Lot Diversity", ylim=c(-18,16))
stripchart(ses ~ diversity, data = SES_ALLdivvl, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(14, 7.4, "A", pos = 4, font = 2, cex = 2.6)


boxplot(ses ~ diversity, data = SES_ALLdivup, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "2019 Urban Prairie Diversity", ylim=c(-18,16))
stripchart(ses ~ diversity, data = SES_ALLdivup, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(14, 7.4, "B", pos = 4, font = 2, cex = 2.6)

dev.off()


library(ggplot2)
library(cowplot)


###Figure 11----

#Length and Origin
SES_lengthandorigin.vl19<- data.frame(
  trmt=rep(SES$trmt,3),
  variable=rep(c("Alien", "Native", "Body Length"), each = nrow(SES)),
  value = c(SES$SES_ori_1, SES$SES_ori_0, SES$SES_bl)
)
SES_lengthandorigin.vl19$variable<-factor(SES_lengthandorigin.vl19$variable, c("Alien", "Native", "Body Length"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots



fig11a<- ggplot(SES_lengthandorigin.vl19, aes(x = variable, y = value, fill = trmt)) +
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

fig11a


#Nesting traits

SES_nesting.vl19<- data.frame(
  trmt=rep(SES$trmt,5),
  variable=rep(c("Wood","Pithy Stems", "Colony", "Cavity", "Soil"), each = nrow(SES)),
  value = c(SES$SES_nest_5, SES$SES_nest_4, SES$SES_nest_3, SES$SES_nest_2, SES$SES_nest_1)
)
SES_nesting.vl19$variable<-factor(SES_nesting.vl19$variable, c("Wood", "Pithy Stems","Colony", "Cavity", "Soil"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots


fig11b<- ggplot(SES_nesting.vl19, aes(x = variable, y = value, fill = trmt)) +
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
fig11b


#Lecty
SES_lecty.vl19<- data.frame(
  trmt=rep(SES$trmt,3),
  variable=rep(c("Specialist","Generalist", "Kleptoparasitic"), each = nrow(SES)),
  value = c(SES$SES_lec_2, SES$SES_lec_1, SES$SES_lec_0)
)
SES_lecty.vl19$variable<-factor(SES_lecty.vl19$variable, c("Specialist", "Generalist","Kleptoparasitic"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots


fig11c<- ggplot(SES_lecty.vl19, aes(x = variable, y = value, fill = trmt)) +
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
fig11c


#Sociality
SES_sociality.vl19<- data.frame(
  trmt=rep(SES$trmt,4),
  variable=rep(c("Parasitic", "Eusocial","Subsocial", "Solitary"), each = nrow(SES)),
  value = c(SES$SES_soc_4, SES$SES_soc_3, SES$SES_soc_2, SES$SES_soc_1)
) #Note that we have to put our funct traits in the reverse order of how we want them to appear on our graph
SES_sociality.vl19$variable<-factor(SES_sociality.vl19$variable, c("Parasitic", "Eusocial","Subsocial", "Solitary"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots

fig11d<- ggplot(SES_sociality.vl19, aes(x = variable, y = value, fill = trmt)) +
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
fig11d


fig11 <- plot_grid(fig11a,fig11b,fig11c,fig11d, labels= c('A','B','C','D'), label_size= 30)
png("Figures/Figure 11 all MP local sp traits.png", width = 1500, height = 1000, pointsize = 20)

fig11

dev.off()



###Figure 12----
treatmentcleanup <- c(VL="Vacant lots", UP="Pocket prairies")
SES$treatmentspelledout <- as.character(treatmentcleanup[SES$trmt])


png("Figures/Figure 12.png", width = 1500, height = 1000, pointsize = 20)

par(mfrow=c(1,2)) # indicates one row, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure



# Parasitic
boxplot(SES_soc_4 ~ treatmentspelledout, data = SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-1.84,1.94), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Parasitic")
stripchart(SES_soc_4 ~ trmt, data = SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.7, 2.42, "A", pos = 4, font = 2, cex = 2)

# Soil Nesting
boxplot(SES_nest_1 ~ treatmentspelledout, data = SES, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-1.84,1.94), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Soil Nesting")
stripchart(SES_nest_1 ~ trmt, data = SES, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.7, 2.42, "B", pos = 4, font = 2, cex = 2)

# Kleptoparasitic
#boxplot(SES_lec_0 ~ treatmentspelledout, data = SES, col = viridis(3, alpha = 0.6),
        #xlab = "Standardized Effect Sizes (SES)", ylab = "", 
       # ylim = c(-2,2), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
       # horizontal = TRUE, las = 1, range = 0, main = "Kleptoparasitic")
#stripchart(SES_lec_0 ~ trmt, data = SES, col = viridis(3),
        #   pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
#abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.87, 2.42, "C", pos = 4, font = 2, cex = 2)

dev.off()
