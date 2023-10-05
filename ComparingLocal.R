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
#Edited again C Shepard 27 Sept 2023
###################################################################################
#Creating the datasets----

t <- read.csv("./btraits_23.csv", row.names=1)
a <- read.csv("./bcomm_23.localanalysis.csv", row.names=1)


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
#SEPT EDIT:Creating two new datasets. One with Farm and Control, one with T1 and T8
a.FS <- rbind(farm, Control)
a.KT<- rbind(T1, T8)
str(a.FS)
str(a.KT)

a.FS <- a.FS[2:362]
a.KT <- a.KT[2:362]
str(a.FS)
str(a.KT)
rowSums(a1[2:362])

# create a vector with the column sums for each species
# species not found in Cleveland will have a 0
sp.FS <- colSums(a.FS)
sp.KT <- colSums(a.KT)
sp.FS
sp.KT


str(t)
names(t)
colnames(t) <- c("bl", "lec", "nest", "soc", "ori")
names(t)

# add the sp vector as a column in the trait matrix, shows which species
# are found in Cleveland and which are absent (i.e. with a 0)
t.FS<-t
t.FS$sp <- sp.FS
t.FS

t.KT<-t
t.KT$sp <- sp.KT
t.KT

# removes any columns (i.e. species) that are not found in Cleveland
a.FS <- a.FS[, colSums(a.FS != 0) > 0]
colSums(a.FS)

a.KT <- a.KT[, colSums(a.KT != 0) > 0]
colSums(a.KT)

# uses the sp values to remove rows of species not collected in Cleveland
# then remove the column because we don't need it anymore
t.FS <- t.FS[t.FS$sp != 0, ]
t.FS <- t.FS[,-6]

t.KT <- t.KT[t.KT$sp != 0, ]
t.KT <- t.KT[,-6]

plot(t.FS)
cor(t.FS, method = c("pearson"), use = "complete.obs")
str(t.FS)

t1.FS <- t.FS #save the original dataset

t.FS$lec <- as.factor(t.FS$lec)
t.FS$nest <- as.factor(t.FS$nest)
t.FS$soc <- as.factor(t.FS$soc)

str(t.FS) # have to keep origin as a integer for the trait distance matrix to work

#check body length for normality
hist(t.FS$bl)
hist(log(t.FS$bl))

t2.FS <- t.FS #create another duplicate dataset before we transform
t2.FS$bl <- log(t2.FS$bl + 1)
str(t2.FS)

#Double check that all species are present in both datasets
#Double check if a species is present in one dataset but not the other
setdiff(colnames(a.FS), rownames(t2.FS))
setdiff(rownames (t2.FS), colnames(a.FS))

rownames(t2.FS) == colnames(a.FS) # we are good to go!

#NOW AGAIN FOR DATA SET 2
plot(t.KT)
cor(t.KT, method = c("pearson"), use = "complete.obs")
str(t.KT)

t1.KT <- t.KT #save the original dataset

t.KT$lec <- as.factor(t.KT$lec)
t.KT$nest <- as.factor(t.KT$nest)
t.KT$soc <- as.factor(t.KT$soc)

str(t.KT) # have to keep origin as a integer for the trait distance matrix to work

#check body length for normality
hist(t.KT$bl)
hist(log(t.KT$bl))

t2.KT <- t.KT #create another duplicate dataset before we transform
t2.KT$bl <- log(t2.KT$bl + 1)
str(t2.KT)

#Double check that all species are present in both datasets
#Double check if a species is present in one dataset but not the other
setdiff(colnames(a.KT), rownames(t2.KT))
setdiff(rownames (t2.KT), colnames(a.KT))

rownames(t2.KT) == colnames(a.KT) # we are good to go!

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

rowSums(a.FS)
rowSums(a.KT)

#observed CWM----
###FS ROUND----
cwm.obs.FS <- functcomp(t2.FS, as.matrix(a.FS), CWM.type = "all")
cwm.obs.FS

#observed taxonomic beta diversity
# create beta part object for analyses
str(a.FS)
b.core.FS <- betapart.core(a.FS)

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
b.dist.FS <- beta.pair(b.core.FS, index.family = "sorensen")
str(b.dist.FS)

#observed functional beta diversity
#create distance matrix with the traits
#optimized feature helps to weight the traits equally
tdis.FS <- gawdis(t2.FS, w.type = "optimized", opti.maxiter = 500)
attr(tdis.FS, "correls")
attr(tdis.FS, "weights")

# save trait weights for the null model
wt.FS <- c(0.38, 0.19, 0.11, 0.12, 0.20)

#now run a principal coordinates analysis (PCoA) so we can collapse these traits into 
#a few continuous axes for the functional diversity calculations
pcoB.FS <- dudi.hillsmith(as.matrix(tdis.FS), scannf = FALSE, nf = 4)
pcoB.FS

# check correlations among axes and traits
cor(pcoB.FS$li, t1.FS, use = "complete.obs")


sum(pcoB.FS$eig[1:4]) / sum(pcoB.FS$eig)
sum(pcoB.FS$eig[1:3]) / sum(pcoB.FS$eig)
sum(pcoB.FS$eig[1:2]) / sum(pcoB.FS$eig)

t.ax.FS <- as.matrix(pcoB.FS$li[1:3])
b.fun.FS <- functional.beta.pair(a.FS, t.ax.FS, index.family = "sorensen")
str(b.fun.FS)

###observed CWM TURO ROUND----
cwm.obs.KT <- functcomp(t2.KT, as.matrix(a.KT), CWM.type = "all")
cwm.obs.KT

#observed taxonomic beta diversity
# create beta part object for analyses
str(a.KT)
b.core.KT <- betapart.core(a.KT)

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
b.dist.KT <- beta.pair(b.core.KT, index.family = "sorensen")
str(b.dist.KT)

#observed functional beta diversity
#create distance matrix with the traits
#optimized feature helps to weight the traits equally
tdis.KT <- gawdis(t2.KT, w.type = "optimized", opti.maxiter = 500)
attr(tdis.KT, "correls")
attr(tdis.KT, "weights")

# save trait weights for the null model
wt.KT <- c(0.38, 0.22, 0.10, 0.11, 0.19)

#now run a principal coordinates analysis (PCoA) so we can collapse these traits into 
#a few continuous axes for the functional diversity calculations
pcoB.KT <- dudi.hillsmith(as.matrix(tdis.KT), scannf = FALSE, nf = 4)
pcoB.KT

# check correlations among axes and traits
cor(pcoB.KT$li, t1.KT, use = "complete.obs")
 


sum(pcoB.KT$eig[1:4]) / sum(pcoB.KT$eig)
sum(pcoB.KT$eig[1:3]) / sum(pcoB.KT$eig)
sum(pcoB.KT$eig[1:2]) / sum(pcoB.KT$eig)

t.ax.KT <- as.matrix(pcoB.KT$li[1:3])
b.fun.KT <- functional.beta.pair(a.KT, t.ax.KT, index.family = "sorensen")
str(b.fun.KT)


#observed functional alpha diversity
#run the rao function first!
bb.rao.FS <- Rao(sample = t(a.FS), dfunc = tdis.FS, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
falpha.FS <- bb.rao.FS$FD$Alpha
falpha.FS


bb.rao.KT <- Rao(sample = t(a.KT), dfunc = tdis.KT, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
falpha.KT <- bb.rao.KT$FD$Alpha
falpha.KT
#We have now calculated all the indices with our observed bee data
#next, we need to run the null modelS!----

#using independent swap method for randomizing the presence/absence matrix
#this will constrain the null communities by species richness and species frequency
#run the null models with 999 iterations
numberReps <- 999


#create matrices to store the results of each iteration of the null model, for each trait and index:----
##FS----
# for cwms
nbl.FS <- nlec_0.FS <- nlec_1.FS <- nlec_2.FS <- nnest_1.FS <- nnest_2.FS <- nnest_3.FS <- nnest_4.FS <- nnest_5.FS <- nsoc_1.FS <- nsoc_2.FS <- nsoc_3.FS <- nsoc_4.FS <- nori_0.FS <- nori_1.FS <- matrix(NA,
                                                                                                                                                               nrow = nrow(a.FS), ncol = numberReps, dimnames = list(rownames(a.FS), paste0("n", 1:numberReps)))

# for taxonomic beta diversity
nbsim.FS <- nbsne.FS <- nbsor.FS <- matrix(NA, nrow = nrow(a.FS), ncol = numberReps, 
                                  dimnames = list(rownames(a.FS), paste0("n", 1:numberReps)))

# for functional beta diversity
nfalpha.FS <- nfsim.FS <- nfsne.FS <- nfsor.FS <- matrix(NA, nrow = nrow(a.FS), ncol = numberReps, 
                                  dimnames = list(rownames(a.FS), paste0("n", 1:numberReps)))

##KT----
# for cwms
nbl.KT <- nlec_0.KT <- nlec_1.KT <- nlec_2.KT <- nnest_1.KT <- nnest_2.KT <- nnest_3.KT <- nnest_4.KT <- nnest_5.KT <- nsoc_1.KT <- nsoc_2.KT <- nsoc_3.KT <- nsoc_4.KT <- nori_0.KT <- nori_1.KT <- matrix(NA,
                                                                                                                                                               nrow = nrow(a.KT), ncol = numberReps, dimnames = list(rownames(a.KT), paste0("n", 1:numberReps)))

# for taxonomic beta diversity
nbsim.KT <- nbsne.KT <- nbsor.KT <- matrix(NA, nrow = nrow(a.KT), ncol = numberReps, 
                                  dimnames = list(rownames(a.KT), paste0("n", 1:numberReps)))

# for functional beta diversity
nfalpha.KT <- nfsim.KT <- nfsne.KT <- nfsor.KT <- matrix(NA, nrow = nrow(a.KT), ncol = numberReps, 
                                             dimnames = list(rownames(a.KT), paste0("n", 1:numberReps)))

#create null model for each repetition:----
##FS EDITION-----

for(i in 1:numberReps){
  print(i) 
  
  # randomized trait matrix
  ntraits.FS <- t2.FS[sample(1:nrow(t2.FS)),]
  rownames(ntraits.FS) <- rownames(t2.FS)
  
  # randomized presence/absence matrix
  nsp.FS <- randomizeMatrix(samp = a.FS, null.model = "independentswap")
  
  # randomized trait distance matrix
  ntdis.FS <- gawdis(ntraits.FS, w.type = "user", W = wt.FS)
  
  # CWM calculations
  cwm.null.FS <- functcomp(x = ntraits.FS, a = as.matrix(nsp.FS), CWM.type = "all")
  nbl.FS[,i] <- cwm.null.FS$bl
  nlec_0.FS[,i] <- cwm.null.FS$lec_0
  nlec_1.FS[,i] <- cwm.null.FS$lec_1
  nlec_2.FS[,i] <- cwm.null.FS$lec_2
  nnest_1.FS[,i] <- cwm.null.FS$nest_1
  nnest_2.FS[,i] <- cwm.null.FS$nest_2
  nnest_3.FS[,i] <- cwm.null.FS$nest_3
  nnest_4.FS[,i] <- cwm.null.FS$nest_4
  nnest_5.FS[,i] <- cwm.null.FS$nest_5
  nsoc_1.FS[,i] <- cwm.null.FS$soc_1
  nsoc_2.FS[,i] <- cwm.null.FS$soc_2
  nsoc_3.FS[,i] <- cwm.null.FS$soc_3
  nsoc_4.FS[,i] <- cwm.null.FS$soc_4
  nori_0.FS[,i] <- cwm.null.FS$ori_0
  nori_1.FS[,i] <- cwm.null.FS$ori_1
  
  # Functional alpha diversity
  nrao.FS <- Rao(sample = t(nsp.FS), dfunc = ntdis.FS, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
  nfalpha.FS[,i] <- nrao.FS$FD$Alpha
  
  # Taxonomic beta diversity indices
  nb.core.FS <- betapart.core(nsp.FS)
  nb.dist.FS <- beta.pair(nb.core.FS, index.family = "sorensen")
  nsim.dist.FS <- as.matrix(nb.dist.FS$beta.sim)
  nsne.dist.FS <- as.matrix(nb.dist.FS$beta.sne)
  nsor.dist.FS <- as.matrix(nb.dist.FS$beta.sor)
  nbsim.FS[,i] <- colMeans(nsim.dist.FS)
  nbsne.FS[,i] <- colMeans(nsne.dist.FS)
  nbsor.FS[,i] <- colMeans(nsor.dist.FS)
  
  # Functional beta diversity indices
  npco.FS <- dudi.hillsmith(as.matrix(ntdis.FS), scannf = FALSE, nf = 3)
  nt.FS <- as.matrix(npco.FS$li)
  nb.fun.FS <- functional.beta.pair(nsp.FS, nt.FS, index.family = "sorensen")
  nfsim.dist.FS <- as.matrix(nb.fun.FS$funct.beta.sim)
  nfsne.dist.FS <- as.matrix(nb.fun.FS$funct.beta.sne)
  nfsor.dist.FS <- as.matrix(nb.fun.FS$funct.beta.sor)
  nfsim.FS[,i] <- colMeans(nfsim.dist.FS)
  nfsne.FS[,i] <- colMeans(nfsne.dist.FS)
  nfsor.FS[,i] <- colMeans(nfsor.dist.FS)
  
}

# save the output matrices
write.csv(nbl.FS, file = "nbl.FS.csv")
write.csv(nlec_0.FS, file = "nlec_0.FS.csv")
write.csv(nlec_1.FS, file = "nlec_1.FS.csv")
write.csv(nlec_2.FS, file = "nlec_2.FS.csv")
write.csv(nnest_1.FS, file = "nnest_1.FS.csv")
write.csv(nnest_2.FS, file = "nnest_2.FS.csv")
write.csv(nnest_3.FS, file = "nnest_3.FS.csv")
write.csv(nnest_4.FS, file = "nnest_4.FS.csv")
write.csv(nnest_5.FS, file = "nnest_5.FS.csv")
write.csv(nsoc_1.FS, file = "nsoc_1.FS.csv")
write.csv(nsoc_2.FS, file = "nsoc_2.FS.csv")
write.csv(nsoc_3.FS, file = "nsoc_3.FS.csv")
write.csv(nsoc_4.FS, file = "nsoc_4.FS.csv")
write.csv(nori_0.FS, file = "nori_0.FS.csv")
write.csv(nori_1.FS, file = "nori_1.FS.csv")

write.csv(nbsim.FS, file = "tbeta_sim.FS.csv")
write.csv(nbsne.FS, file = "tbeta_sne.FS.csv")
write.csv(nbsor.FS, file = "tbeta_sor.FS.csv")

write.csv(nfalpha.FS, file = "falpha.FS.csv")

write.csv(nfsim.FS, file = "fbeta_sim.FS.csv")
write.csv(nfsne.FS, file = "fbeta_sne.FS.csv")
write.csv(nfsor.FS, file = "fbeta_sor.FS.csv")

# load the output matrices

#I manually moved all of the csv files to the folder Urban to Local_Null Models_HillSmith_3 Axes.FS
nbl.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nbl.FS.csv", row.names=1)
nlec_0.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nlec_0.FS.csv", row.names=1)
nlec_1.FS<- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nlec_1.FS.csv", row.names=1)
nlec_2.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nlec_2.FS.csv", row.names=1)
nnest_1.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nnest_1.FS.csv", row.names=1)
nnest_2.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nnest_2.FS.csv", row.names=1)
nnest_3.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nnest_3.FS.csv", row.names=1)
nnest_4.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nnest_4.FS.csv", row.names=1)
nnest_5.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nnest_5.FS.csv", row.names=1)
nsoc_1.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nsoc_1.FS.csv", row.names=1)
nsoc_2.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nsoc_2.FS.csv", row.names=1)
nsoc_3.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nsoc_3.FS.csv", row.names=1)
nsoc_4.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nsoc_4.FS.csv", row.names=1)
nori_0.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nori_0.FS.csv", row.names=1)
nori_1.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/nori_1.FS.csv", row.names=1)

nfalpha.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/falpha.FS.csv", row.names=1)

nbsim.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/tbeta_sim.FS.csv", row.names=1)
nbsne.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/tbeta_sne.FS.csv", row.names=1)
nbsor.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/tbeta_sim.FS.csv", row.names=1)

nfsim.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/fbeta_sim.FS.csv", row.names=1)
nfsne.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/fbeta_sne.FS.csv", row.names=1)
nfsor.FS <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.FS/fbeta_sor.FS.csv", row.names=1)


##KT EDITION-----

for(i in 1:numberReps){
  print(i) 
  
  # randomized trait matrix
  ntraits.KT <- t2.KT[sample(1:nrow(t2.KT)),]
  rownames(ntraits.KT) <- rownames(t2.KT)
  
  # randomized presence/absence matrix
  nsp.KT <- randomizeMatrix(samp = a.KT, null.model = "independentswap")
  
  # randomized trait distance matrix
  ntdis.KT <- gawdis(ntraits.KT, w.type = "user", W = wt.KT)
  
  # CWM calculations
  cwm.null.KT <- functcomp(x = ntraits.KT, a = as.matrix(nsp.KT), CWM.type = "all")
  nbl.KT[,i] <- cwm.null.KT$bl
  nlec_0.KT[,i] <- cwm.null.KT$lec_0
  nlec_1.KT[,i] <- cwm.null.KT$lec_1
  nlec_2.KT[,i] <- cwm.null.KT$lec_2
  nnest_1.KT[,i] <- cwm.null.KT$nest_1
  nnest_2.KT[,i] <- cwm.null.KT$nest_2
  nnest_3.KT[,i] <- cwm.null.KT$nest_3
  nnest_4.KT[,i] <- cwm.null.KT$nest_4
  nnest_5.KT[,i] <- cwm.null.KT$nest_5
  nsoc_1.KT[,i] <- cwm.null.KT$soc_1
  nsoc_2.KT[,i] <- cwm.null.KT$soc_2
  nsoc_3.KT[,i] <- cwm.null.KT$soc_3
  nsoc_4.KT[,i] <- cwm.null.KT$soc_4
  nori_0.KT[,i] <- cwm.null.KT$ori_0
  nori_1.KT[,i] <- cwm.null.KT$ori_1
  
  # Functional alpha diversity
  nrao.KT <- Rao(sample = t(nsp.KT), dfunc = ntdis.KT, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL)
  nfalpha.KT[,i] <- nrao.KT$FD$Alpha
  
  # Taxonomic beta diversity indices
  nb.core.KT <- betapart.core(nsp.KT)
  nb.dist.KT <- beta.pair(nb.core.KT, index.family = "sorensen")
  nsim.dist.KT <- as.matrix(nb.dist.KT$beta.sim)
  nsne.dist.KT <- as.matrix(nb.dist.KT$beta.sne)
  nsor.dist.KT <- as.matrix(nb.dist.KT$beta.sor)
  nbsim.KT[,i] <- colMeans(nsim.dist.KT)
  nbsne.KT[,i] <- colMeans(nsne.dist.KT)
  nbsor.KT[,i] <- colMeans(nsor.dist.KT)
  
  # Functional beta diversity indices
  npco.KT <- dudi.hillsmith(as.matrix(ntdis.KT), scannf = FALSE, nf = 3)
  nt.KT <- as.matrix(npco.KT$li)
  nb.fun.KT <- functional.beta.pair(nsp.KT, nt.KT, index.family = "sorensen")
  nfsim.dist.KT <- as.matrix(nb.fun.KT$funct.beta.sim)
  nfsne.dist.KT <- as.matrix(nb.fun.KT$funct.beta.sne)
  nfsor.dist.KT <- as.matrix(nb.fun.KT$funct.beta.sor)
  nfsim.KT[,i] <- colMeans(nfsim.dist.KT)
  nfsne.KT[,i] <- colMeans(nfsne.dist.KT)
  nfsor.KT[,i] <- colMeans(nfsor.dist.KT)
  
}

# save the output matrices
write.csv(nbl.KT, file = "nbl.KT.csv")
write.csv(nlec_0.KT, file = "nlec_0.KT.csv")
write.csv(nlec_1.KT, file = "nlec_1.KT.csv")
write.csv(nlec_2.KT, file = "nlec_2.KT.csv")
write.csv(nnest_1.KT, file = "nnest_1.KT.csv")
write.csv(nnest_2.KT, file = "nnest_2.KT.csv")
write.csv(nnest_3.KT, file = "nnest_3.KT.csv")
write.csv(nnest_4.KT, file = "nnest_4.KT.csv")
write.csv(nnest_5.KT, file = "nnest_5.KT.csv")
write.csv(nsoc_1.KT, file = "nsoc_1.KT.csv")
write.csv(nsoc_2.KT, file = "nsoc_2.KT.csv")
write.csv(nsoc_3.KT, file = "nsoc_3.KT.csv")
write.csv(nsoc_4.KT, file = "nsoc_4.KT.csv")
write.csv(nori_0.KT, file = "nori_0.KT.csv")
write.csv(nori_1.KT, file = "nori_1.KT.csv")

write.csv(nbsim.KT, file = "tbeta_sim.KT.csv")
write.csv(nbsne.KT, file = "tbeta_sne.KT.csv")
write.csv(nbsor.KT, file = "tbeta_sor.KT.csv")

write.csv(nfalpha.KT, file = "falpha.KT.csv")

write.csv(nfsim.KT, file = "fbeta_sim.KT.csv")
write.csv(nfsne.KT, file = "fbeta_sne.KT.csv")
write.csv(nfsor.KT, file = "fbeta_sor.KT.csv")

# load the output matrices

#I manually moved all of the csv files to the folder Urban to Local_Null Models_HillSmith_3 Axes.KT
nbl.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nbl.KT.csv", row.names=1)
nlec_0.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nlec_0.KT.csv", row.names=1)
nlec_1.KT<- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nlec_1.KT.csv", row.names=1)
nlec_2.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nlec_2.KT.csv", row.names=1)
nnest_1.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nnest_1.KT.csv", row.names=1)
nnest_2.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nnest_2.KT.csv", row.names=1)
nnest_3.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nnest_3.KT.csv", row.names=1)
nnest_4.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nnest_4.KT.csv", row.names=1)
nnest_5.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nnest_5.KT.csv", row.names=1)
nsoc_1.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nsoc_1.KT.csv", row.names=1)
nsoc_2.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nsoc_2.KT.csv", row.names=1)
nsoc_3.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nsoc_3.KT.csv", row.names=1)
nsoc_4.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nsoc_4.KT.csv", row.names=1)
nori_0.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nori_0.KT.csv", row.names=1)
nori_1.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/nori_1.KT.csv", row.names=1)

nfalpha.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/falpha.KT.csv", row.names=1)

nbsim.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/tbeta_sim.KT.csv", row.names=1)
nbsne.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/tbeta_sne.KT.csv", row.names=1)
nbsor.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/tbeta_sim.KT.csv", row.names=1)

nfsim.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/fbeta_sim.KT.csv", row.names=1)
nfsne.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/fbeta_sne.KT.csv", row.names=1)
nfsor.KT <- read.csv("Urban to Local_Null Models_HillSmith_3 Axes.KT/fbeta_sor.KT.csv", row.names=1)


# SES Calculations----
#calculate standardized effect sizes (SES) for each trait and index
#the effect size is the difference between the observed value and the expected one
#then divide the effect size by the standard deviation of the null distribution to get the standardized effect size
#allows comparison among sites with different numbers of species


## community weighted means----
###FS VERSION----
## body length
SES_bl.FS <- (cwm.obs.FS$bl - apply(nbl.FS, MARGIN = 1, mean)) / apply(nbl.FS, MARGIN = 1, sd, na.rm=T)
SES_bl.FS

## lec_0 - Kleptoparasitic
SES_lec_0.FS <- (cwm.obs.FS$lec_0 - apply(nlec_0.FS, MARGIN = 1, mean)) / apply(nlec_0.FS, MARGIN = 1, sd, na.rm=T)
SES_lec_0.FS

## lec_1 - Generalist
SES_lec_1.FS <- (cwm.obs.FS$lec_1 - apply(nlec_1.FS, MARGIN = 1, mean)) / apply(nlec_1.FS, MARGIN = 1, sd, na.rm=T)
SES_lec_1.FS

## lec_2 - Specialist
SES_lec_2.FS <- (cwm.obs.FS$lec_2 - apply(nlec_2.FS, MARGIN = 1, mean)) / apply(nlec_2.FS, MARGIN = 1, sd, na.rm=T)
SES_lec_2.FS

## nest_1 - Soil
SES_nest_1.FS <- (cwm.obs.FS$nest_1 - apply(nnest_1.FS, MARGIN = 1, mean)) / apply(nnest_1.FS, MARGIN = 1, sd, na.rm=T)
SES_nest_1.FS

## nest_2 - Cavity
SES_nest_2.FS<- (cwm.obs.FS$nest_2 - apply(nnest_2.FS, MARGIN = 1, mean)) / apply(nnest_2.FS, MARGIN = 1, sd, na.rm=T)
SES_nest_2.FS

## nest_3 - Hive
SES_nest_3.FS <- (cwm.obs.FS$nest_3 - apply(nnest_3.FS, MARGIN = 1, mean)) / apply(nnest_3.FS, MARGIN = 1, sd, na.rm=T)
SES_nest_3.FS

## nest_4 - Pithy Stems
SES_nest_4.FS <- (cwm.obs.FS$nest_4 - apply(nnest_4.FS, MARGIN = 1, mean)) / apply(nnest_4.FS, MARGIN = 1, sd, na.rm=T)
SES_nest_4.FS

## nest_5 - Wood
SES_nest_5.FS <- (cwm.obs.FS$nest_5 - apply(nnest_5.FS, MARGIN = 1, mean)) / apply(nnest_5.FS, MARGIN = 1, sd, na.rm=T)
SES_nest_5.FS

## soc_1 - Subsocial
SES_soc_1.FS <- (cwm.obs.FS$soc_1 - apply(nsoc_1, MARGIN = 1, mean)) / apply(nsoc_1.FS, MARGIN = 1, sd, na.rm=T)
SES_soc_1.FS

## soc_2 - Solitary
SES_soc_2.FS <- (cwm.obs.FS$soc_2 - apply(nsoc_2.FS, MARGIN = 1, mean)) / apply(nsoc_2.FS, MARGIN = 1, sd, na.rm=T)
SES_soc_2.FS

## soc_3 - Eusocial
SES_soc_3.FS <- (cwm.obs.FS$soc_3 - apply(nsoc_3.FS, MARGIN = 1, mean)) / apply(nsoc_3.FS, MARGIN = 1, sd, na.rm=T)
SES_soc_3.FS

## soc_4 - Parasitic
SES_soc_4.FS <- (cwm.obs.FS$soc_4 - apply(nsoc_4.FS, MARGIN = 1, mean)) / apply(nsoc_4.FS, MARGIN = 1, sd, na.rm=T)
SES_soc_4.FS

## ori_0 - Native
SES_ori_0.FS <- (cwm.obs.FS$ori_0 - apply(nori_0.FS, MARGIN = 1, mean)) / apply(nori_0.FS, MARGIN = 1, sd, na.rm=T)
SES_ori_0.FS

## ori_2 - Exotic
SES_ori_1.FS <- (cwm.obs.FS$ori_1 - apply(nori_1.FS, MARGIN = 1, mean)) / apply(nori_1.FS, MARGIN = 1, sd, na.rm=T)
SES_ori_1.FS

## taxonomic beta diversity
beta.sor.FS <- as.matrix(b.dist.FS$beta.sor)
beta.sor.FS <- colMeans(beta.sor.FS)

beta.sim.FS <- as.matrix(b.dist.FS$beta.sim)
beta.sim.FS <- colMeans(beta.sim.FS)

beta.sne.FS <- as.matrix(b.dist.FS$beta.sne)
beta.sne.FS <- colMeans(beta.sne.FS)

beta.t.FS <- data.frame(beta.sor.FS, beta.sim.FS, beta.sne.FS)

## taxonomic diveristy - beta sor
SES_bsor.FS <- (beta.t.FS$beta.sor - apply(nbsor.FS, MARGIN = 1, mean)) / apply(nbsor.FS, MARGIN = 1, sd, na.rm=T)
SES_bsor.FS

## taxonomic diveristy - beta sim
SES_bsim.FS <- (beta.t.FS$beta.sim - apply(nbsim.FS, MARGIN = 1, mean)) / apply(nbsim.FS, MARGIN = 1, sd, na.rm=T)
SES_bsim.FS

## taxonomic diveristy - beta sne
SES_bsne.FS <- (beta.t.FS$beta.sne - apply(nbsne.FS, MARGIN = 1, mean)) / apply(nbsne.FS, MARGIN = 1, sd, na.rm=T)
SES_bsne.FS

## functional diversity
fbeta.sor.FS <- as.matrix(b.fun.FS$funct.beta.sor)
fbeta.sor.FS <- colMeans(fbeta.sor.FS)

fbeta.sim.FS <- as.matrix(b.fun.FS$funct.beta.sim)
fbeta.sim.FS <- colMeans(fbeta.sim.FS)

fbeta.sne.FS <- as.matrix(b.fun.FS$funct.beta.sne)
fbeta.sne.FS <- colMeans(fbeta.sne.FS)

beta.f.FS <- data.frame(falpha.FS, fbeta.sor.FS, fbeta.sim.FS, fbeta.sne.FS)

## functional alpha  - Rao
SES_falpha.FS <- (beta.f.FS$falpha - apply(nfalpha.FS, MARGIN = 1, mean)) / apply(nfalpha.FS, MARGIN = 1, sd, na.rm=T)
SES_falpha.FS

## functional diveristy - beta sor
SES_fbsor.FS <- (beta.f.FS$fbeta.sor - apply(nfsor.FS, MARGIN = 1, mean)) / apply(nfsor.FS, MARGIN = 1, sd, na.rm=T)
SES_fbsor.FS

## functional diveristy - beta sim
SES_fbsim.FS <- (beta.f.FS$fbeta.sim - apply(nfsim.FS, MARGIN = 1, mean)) / apply(nfsim.FS, MARGIN = 1, sd, na.rm=T)
SES_fbsim.FS

## functional diveristy - beta sne
SES_fbsne.FS <- (beta.f.FS$fbeta.sne - apply(nfsne.FS, MARGIN = 1, mean)) / apply(nfsne.FS, MARGIN = 1, sd, na.rm=T)
SES_fbsne.FS

## combine all indices into one matrix
SES.all.FS <- as.data.frame(cbind(SES_bl.FS, SES_lec_0.FS, SES_lec_1.FS, SES_lec_2.FS, SES_ori_0.FS, SES_ori_1.FS, 
                           SES_nest_1.FS, SES_nest_2.FS, SES_nest_3.FS, SES_nest_4.FS, SES_nest_5.FS,
                           SES_soc_1.FS, SES_soc_2.FS, SES_soc_3.FS, SES_soc_4.FS, SES_bsor.FS, SES_bsim.FS,
                           SES_bsne.FS, SES_falpha.FS, SES_fbsor.FS, SES_fbsim.FS, SES_fbsne.FS))

write.csv(SES.all.FS, file = "SES_Local.FS.csv")
#import the SES data
SES.all.FS <- read.csv("SES_Local.FS.csv", row.names = 1)

###KT VERSION----
## body length
SES_bl.KT <- (cwm.obs.KT$bl - apply(nbl.KT, MARGIN = 1, mean)) / apply(nbl.KT, MARGIN = 1, sd, na.rm=T)
SES_bl.KT

## lec_0 - Kleptoparasitic
SES_lec_0.KT <- (cwm.obs.KT$lec_0 - apply(nlec_0.KT, MARGIN = 1, mean)) / apply(nlec_0.KT, MARGIN = 1, sd, na.rm=T)
SES_lec_0.KT

## lec_1 - Generalist
SES_lec_1.KT <- (cwm.obs.KT$lec_1 - apply(nlec_1.KT, MARGIN = 1, mean)) / apply(nlec_1.KT, MARGIN = 1, sd, na.rm=T)
SES_lec_1.KT

## lec_2 - Specialist
SES_lec_2.KT <- (cwm.obs.KT$lec_2 - apply(nlec_2.KT, MARGIN = 1, mean)) / apply(nlec_2.KT, MARGIN = 1, sd, na.rm=T)
SES_lec_2.KT

## nest_1 - Soil
SES_nest_1.KT <- (cwm.obs.KT$nest_1 - apply(nnest_1.KT, MARGIN = 1, mean)) / apply(nnest_1.KT, MARGIN = 1, sd, na.rm=T)
SES_nest_1.KT

## nest_2 - Cavity
SES_nest_2.KT <- (cwm.obs.KT$nest_2 - apply(nnest_2.KT, MARGIN = 1, mean)) / apply(nnest_2.KT, MARGIN = 1, sd, na.rm=T)
SES_nest_2.KT

## nest_3 - Hive
SES_nest_3.KT <- (cwm.obs.KT$nest_3 - apply(nnest_3.KT, MARGIN = 1, mean)) / apply(nnest_3.KT, MARGIN = 1, sd, na.rm=T)
SES_nest_3.KT

## nest_4 - Pithy Stems
SES_nest_4.KT <- (cwm.obs.KT$nest_4 - apply(nnest_4.KT, MARGIN = 1, mean)) / apply(nnest_4.KT, MARGIN = 1, sd, na.rm=T)
SES_nest_4.KT

## nest_5 - Wood
SES_nest_5.KT <- (cwm.obs.KT$nest_5 - apply(nnest_5.KT, MARGIN = 1, mean)) / apply(nnest_5.KT, MARGIN = 1, sd, na.rm=T)
SES_nest_5.KT

## soc_1 - Subsocial
SES_soc_1.KT <- (cwm.obs.KT$soc_1 - apply(nsoc_1.KT, MARGIN = 1, mean)) / apply(nsoc_1.KT, MARGIN = 1, sd, na.rm=T)
SES_soc_1.KT

## soc_2 - Solitary
SES_soc_2.KT <- (cwm.obs.KT$soc_2 - apply(nsoc_2.KT, MARGIN = 1, mean)) / apply(nsoc_2.KT, MARGIN = 1, sd, na.rm=T)
SES_soc_2.KT

## soc_3 - Eusocial
SES_soc_3.KT <- (cwm.obs.KT$soc_3 - apply(nsoc_3.KT, MARGIN = 1, mean)) / apply(nsoc_3.KT, MARGIN = 1, sd, na.rm=T)
SES_soc_3.KT

## soc_4 - Parasitic
SES_soc_4.KT <- (cwm.ob.KTs$soc_4 - apply(nsoc_4.KT, MARGIN = 1, mean)) / apply(nsoc_4.KT, MARGIN = 1, sd, na.rm=T)
SES_soc_4.KT

## ori_0 - Native
SES_ori_0.KT <- (cwm.obs.KT$ori_0 - apply(nori_0.KT, MARGIN = 1, mean)) / apply(nori_0.KT, MARGIN = 1, sd, na.rm=T)
SES_ori_0.KT

## ori_2 - Exotic
SES_ori_1.KT <- (cwm.obs.KT$ori_1 - apply(nori_1.KT, MARGIN = 1, mean)) / apply(nori_1.KT, MARGIN = 1, sd, na.rm=T)
SES_ori_1.KT

## taxonomic beta diversity
beta.sor.KT <- as.matrix(b.dist.KT$beta.sor)
beta.sor.KT <- colMeans(beta.sor.KT)

beta.sim.KT <- as.matrix(b.dist.KT$beta.sim)
beta.sim.KT <- colMeans(beta.sim.KT)

beta.sne.KT <- as.matrix(b.dist.KT$beta.sne)
beta.sne.KT <- colMeans(beta.sne.KT)

beta.t.KT <- data.frame(beta.sor.KT, beta.sim.KT, beta.sne.KT)

## taxonomic diveristy - beta sor
SES_bsor.KT <- (beta.t.KT$beta.sor - apply(nbsor.KT, MARGIN = 1, mean)) / apply(nbsor.KT, MARGIN = 1, sd, na.rm=T)
SES_bsor.KT

## taxonomic diveristy - beta sim
SES_bsim.KT<- (beta.t.KT$beta.sim - apply(nbsim.KT, MARGIN = 1, mean)) / apply(nbsim.KT, MARGIN = 1, sd, na.rm=T)
SES_bsim.KT

## taxonomic diveristy - beta sne
SES_bsne.KT <- (beta.t.KT$beta.sne - apply(nbsne.KT, MARGIN = 1, mean)) / apply(nbsne.KT, MARGIN = 1, sd, na.rm=T)
SES_bsne.KT

## functional diversity
fbeta.sor.KT <- as.matrix(b.fun.KT$funct.beta.sor)
fbeta.sor.KT <- colMeans(fbeta.sor.KT)

fbeta.sim.KT <- as.matrix(b.fun.KT$funct.beta.sim)
fbeta.sim.KT <- colMeans(fbeta.sim.KT)

fbeta.sne.KT <- as.matrix(b.fun.KT$funct.beta.sne)
fbeta.sne.KT <- colMeans(fbeta.sne.KT)

beta.f.KT <- data.frame(falpha.KT, fbeta.sor.KT, fbeta.sim.KT, fbeta.sne.KT)

## functional alpha  - Rao
SES_falpha.KT <- (beta.f.KT$falpha - apply(nfalpha.KT, MARGIN = 1, mean)) / apply(nfalpha.KT, MARGIN = 1, sd, na.rm=T)
SES_falpha.KT

## functional diveristy - beta sor
SES_fbsor.KT <- (beta.f.KT$fbeta.sor - apply(nfsor.KT, MARGIN = 1, mean)) / apply(nfsor.KT, MARGIN = 1, sd, na.rm=T)
SES_fbsor.KT

## functional diveristy - beta sim
SES_fbsim.KT <- (beta.f.KT$fbeta.sim - apply(nfsim.KT, MARGIN = 1, mean)) / apply(nfsim.KT, MARGIN = 1, sd, na.rm=T)
SES_fbsim.KT

## functional diveristy - beta sne
SES_fbsne.KT <- (beta.f.KT$fbeta.sne - apply(nfsne.KT, MARGIN = 1, mean)) / apply(nfsne.KT, MARGIN = 1, sd, na.rm=T)
SES_fbsne.KT

## combine all indices into one matrix
SES.all.KT <- as.data.frame(cbind(SES_bl.KT, SES_lec_0.KT, SES_lec_1.KT, SES_lec_2.KT, SES_ori_0.KT, SES_ori_1.KT, 
                               SES_nest_1.KT, SES_nest_2.KT, SES_nest_3.KT, SES_nest_4.KT, SES_nest_5.KT,
                               SES_soc_1.KT, SES_soc_2.KT, SES_soc_3.KT, SES_soc_4.KT, SES_bsor.KT, SES_bsim.KT,
                               SES_bsne.KT, SES_falpha.KT, SES_fbsor.KT, SES_fbsim.KT, SES_fbsne.KT))

write.csv(SES.all.KT, file = "SES_Local.KT.csv")
#import the SES data
SES.all.KT <- read.csv("SES_Local.KT.csv", row.names = 1)


# import the landscape data-----
land <- read.csv("landscape.localanalysis.csv", row.names = 1)
str(land)

land$trmt <- as.factor(land$trmt)
str(land)

## merge landscape data with SES data
SES <- merge(SES.all, land, by = c("row.names"))
str(SES)

SES$site <- as.factor(SES$Row.names)
str(SES)

write.csv(SES, file = "SES.csv")
## pull out data for each treatment
farm <- SES[which(SES$trmt == "Farm"),]
str(farm)

control <- SES[which(SES$trmt == "Control"),]
str(control)

T1 <- SES[which(SES$trmt == "T1"),]
str(T1)

T8 <- SES[which(SES$trmt == "T8"),]
str(T8)

##################################################################################################
# Make figures!

library(viridis)

# pull out treatments that you want to include
fc <- rbind(farm, control)
str(fc)
fc <- droplevels(fc) #drops T1 and T8
str(fc)

# changes the treatment names for the figure
levels(fc$trmt)[levels(fc$trmt)=='Farm'] <- 'Urban Farm'
levels(fc$trmt)[levels(fc$trmt)=='Control'] <- 'Vacant Lot'
levels(fc$trmt)

# this will save the figure as a png, it will be good quality
png("Figures/Figure 1 v2.png", width = 1500, height = 1000, pointsize = 20)

par(mfrow=c(2,2)) # indicates two rows, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure

# hive nesting
boxplot(SES_nest_3 ~ trmt, data = fc, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-2,3), 
        cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Hive Nesting")
stripchart(SES_ori_0 ~ trmt, data = fc, col = viridis(3),  ylim = c(-2,3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

# native
boxplot(SES_ori_0 ~ trmt, data = fc, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        #ylim = c(-2,3), 
        cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Native")
stripchart(SES_ori_0 ~ trmt, data = fc, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

# specialists
boxplot(SES_lec_2 ~ trmt, data = fc, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        #ylim = c(-2,3), 
        cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Specialists")
stripchart(SES_ori_0 ~ trmt, data = fc, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

# non-native
boxplot(SES_ori_1 ~ trmt, data = fc, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        #ylim = c(-2,3), 
        cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Non-Native")
stripchart(SES_ori_1 ~ trmt, data = fc, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)


dev.off()

##################################################################################################

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

##Below is code for not comparing all at once----
#Don't do this code# FS_bl <- as.data.frame(cbind(farm$SES_bl, control$SES_bl))
#The above makes a matrix, but it doesn't work because the different treatments have different locations. can't be compared equally
#The below takes more time, but it works okay. I am pulling out data site by site to build the matrix I want
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
FS_comps


#yep, you can also do it this way
f.farm <- SES[which(SES$trmt == "Farm"),]
str(f.farm)

f.control <- SES[which(SES$trmt == "Control"),]
str(f.control)

FS_comps <- as.data.frame(rbind(f.control, f.farm))#here you bind by rows rather than columns (like the code above)

#Okay so the above created a matrix of Francis' results
with(FS_comps, bartlett.test(SES_bl ~ trmt))
#####IT WORKS!!!!!
#Okay so to do
#1- make a matrix of Katie results
KT_comps<-SES[3:33,]
KT_comps<-KT_comps[-3,]
KT_comps<-KT_comps[-3,]
KT_comps<-KT_comps[-3,]
KT_comps<-KT_comps[-3,]
KT_comps<-KT_comps[-3,]
KT_comps<-KT_comps[-3,]
KT_comps<-KT_comps[-3,]
KT_comps<-KT_comps[-3,]
KT_comps<-KT_comps[-3,]
KT_comps<-KT_comps[-3,]
KT_comps<-KT_comps[-7,]
KT_comps<-KT_comps[-9,]
KT_comps<-KT_comps[-13,]
KT_comps<-KT_comps[-13,]
KT_comps<-KT_comps[-13,]
#2- Make a matrix of both vacant lot results
VL_comps<-SES[3:32,]
VL_comps<-VL_comps[-2,]
VL_comps<-VL_comps[-2,]
VL_comps<-VL_comps[-12,]
VL_comps<-VL_comps[-13,]
VL_comps<-VL_comps[-20,]
VL_comps<-VL_comps[-20,]
VL_comps<-VL_comps[-20,]
VL_comps<-VL_comps[-20,]
VL_comps<-VL_comps[-13,]
VL_comps<-VL_comps[-14,]
VL_comps<-VL_comps[-14,]
VL_comps<-VL_comps[-18,]
VL_comps<-VL_comps[-15,]
#rewrite code for comparisons to be about each thing

#back to comparisons of body length----
with(FS_comps, bartlett.test(SES_bl ~ trmt))
with(KT_comps, bartlett.test(SES_bl ~ trmt))
with(VL_comps, bartlett.test(SES_bl ~ trmt))
with(FS_comps, ad.test(SES_bl))
with(KT_comps, ad.test(SES_bl))
with(VL_comps, ad.test(SES_bl))

FS_bl.mod.full <- glm(SES_bl ~ trmt + pland + enn, family = gaussian, data = FS_comps)
summary(FS_bl.mod.full)
step(FS_bl.mod.full)

FS_bl.mod.red <- glm(SES_bl ~ pland + enn, family = gaussian, data = FS_comps)
summary(FS_bl.mod.red)
qqnorm(resid(FS_bl.mod.red))
qqline(resid(FS_bl.mod.red))
plot(simulateResiduals(FS_bl.mod.red))
densityPlot(rstudent(FS_bl.mod.red)) # check density estimate of the distribution of residuals
simout <- simulateResiduals(fittedModel = FS_bl.mod.red, plot = T)
plotResiduals(simout, form = FS_comps$pland)
plotResiduals(simout, form = FS_comps$enn)
outlierTest(FS_bl.mod.red)
influenceIndexPlot(FS_bl.mod.red, vars = c("Cook"), id = list(n = 3))
#It doesn't look like there really is an outlier with Bonferroini


FS_bl.mod.null <- glm(SES_bl ~ 1, family = gaussian, data = FS_comps)

# model comparison techniques
anova(FS_bl.mod.full, FS_bl.mod.red, FS_bl.mod.null, test = "F")
AICctab(FS_bl.mod.full, FS_bl.mod.red, FS_bl.mod.null)

Anova(FS_bl.mod.red)

effect_plot(FS_bl.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = '(PLAND)', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(FS_bl.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = '(ENN)', y.label = 'Standardized Effect Sizes (SES)')



KT_bl.mod.full <- glm(SES_bl ~ trmt + pland + enn, family = gaussian, data = KT_comps)
summary(KT_bl.mod.full)
step(KT_bl.mod.full)

KT_bl.mod.null <- glm(SES_bl ~ 1, family = gaussian, data = KT_comps)
summary(KT_bl.mod.null)
qqnorm(resid(KT_bl.mod.null))
qqline(resid(KT_bl.mod.null))
plot(simulateResiduals(KT_bl.mod.null))
densityPlot(rstudent(KT_bl.mod.null)) # check density estimate of the distribution of residuals
outlierTest(KT_bl.mod.null)
influenceIndexPlot(KT_bl.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(KT_bl.mod.full, KT_bl.mod.null, test = "F")
AICctab(KT_bl.mod.full, KT_bl.mod.null)


VL_bl.mod.full <- glm(SES_bl ~ trmt + pland + enn, family = gaussian, data = VL_comps)
summary(VL_bl.mod.full)
step(VL_bl.mod.full)

VL_bl.mod.null <- glm(SES_bl ~ 1, family = gaussian, data = VL_comps)
summary(VL_bl.mod.null)
qqnorm(resid(VL_bl.mod.null))
qqline(resid(VL_bl.mod.null))
plot(simulateResiduals(VL_bl.mod.null))
densityPlot(rstudent(VL_bl.mod.null)) # check density estimate of the distribution of residuals
outlierTest(VL_bl.mod.null)
influenceIndexPlot(VL_bl.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(VL_bl.mod.full, VL_bl.mod.null, test = "F")
AICctab(VL_bl.mod.full, VL_bl.mod.null)



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
with(FS_comps, bartlett.test(SES_lec_0 ~ trmt))
with(KT_comps, bartlett.test(SES_lec_0 ~ trmt))
with(VL_comps, bartlett.test(SES_lec_0 ~ trmt))
with(FS_comps, ad.test(SES_lec_0))
with(KT_comps, ad.test(SES_lec_0))
with(VL_comps, ad.test(SES_lec_0))


FS_lec_0.mod.full <- glm(SES_lec_0 ~ trmt + pland + enn, family = gaussian, data = FS_comps)
summary(FS_lec_0.mod.full)
step(FS_lec_0.mod.full)

FS_lec_0.mod.null <- glm(SES_lec_0 ~ 1, family = gaussian, data = FS_comps)
summary(FS_lec_0.mod.null)
qqnorm(resid(FS_lec_0.mod.null))
qqline(resid(FS_lec_0.mod.null))
plot(simulateResiduals(FS_lec_0.mod.null))
densityPlot(rstudent(FS_lec_0.mod.null)) # check density estimate of the distribution of residuals
outlierTest(FS_lec_0.mod.null)
influenceIndexPlot(FS_lec_0.mod.null, vars = c("Cook"), id = list(n = 3))
#It looks like the first point might be a significant outlier
FS_lec_0.mod.null2 <- update(FS_lec_0.mod.null, subset = -c(1))
summary(FS_lec_0.mod.null2)
compareCoefs(FS_lec_0.mod.null, FS_lec_0.mod.null2) # compares estimated coefficients and their standard errors
#For some reason I got an error when I tried the above code
#And I can't test the below comparisons, because they don't run with it having 1 fewer data point

# model comparison techniques
anova(FS_lec_0.mod.full, FS_lec_0.mod.null,  test = "F")
AICctab(FS_lec_0.mod.full, FS_lec_0.mod.null)

KT_lec_0.mod.full <- glm(SES_lec_0 ~ trmt + pland + enn, family = gaussian, data = KT_comps)
summary(KT_lec_0.mod.full)
step(KT_lec_0.mod.full)

KT_lec_0.mod.red <- glm(SES_lec_0 ~ enn, family = gaussian, data = KT_comps)
summary(KT_lec_0.mod.red)
qqnorm(resid(KT_lec_0.mod.red))
qqline(resid(KT_lec_0.mod.red))
plot(simulateResiduals(KT_lec_0.mod.red))
#I am adding the below code to look at the graphs individually. the second one has red lines
plotQQunif(KT_lec_0.mod.red)
plotResiduals(KT_lec_0.mod.red)
#It says that there are significant quantile deviations

densityPlot(rstudent(KT_lec_0.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_lec_0.mod.red)
influenceIndexPlot(KT_lec_0.mod.red, vars = c("Cook"), id = list(n = 3))

KT_lec_0.mod.null <- glm(SES_lec_0 ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_lec_0.mod.full, KT_lec_0.mod.red, KT_lec_0.mod.null,  test = "F")
AICctab(KT_lec_0.mod.full, KT_lec_0.mod.red, KT_lec_0.mod.null)

VL_lec_0.mod.full <- glm(SES_lec_0 ~ trmt + pland + enn, family = gaussian, data = VL_comps)
summary(VL_lec_0.mod.full)
step(VL_lec_0.mod.full)

VL_lec_0.mod.red <- glm(SES_lec_0 ~ trmt+pland, family = gaussian, data = VL_comps)
summary(VL_lec_0.mod.red)
qqnorm(resid(VL_lec_0.mod.red))
qqline(resid(VL_lec_0.mod.red))
plot(simulateResiduals(VL_lec_0.mod.red))
#I am adding the below code to look at the graphs individually. the second one has red lines
plotQQunif(VL_lec_0.mod.red)
plotResiduals(VL_lec_0.mod.red)
densityPlot(rstudent(VL_lec_0.mod.red)) # check density estimate of the distribution of residuals
outlierTest(VL_lec_0.mod.red)
influenceIndexPlot(VL_lec_0.mod.red, vars = c("Cook"), id = list(n = 3))

VL_lec_0.mod.null <- glm(SES_lec_0 ~ 1, family = gaussian, data = VL_comps)

# model comparison techniques
anova(VL_lec_0.mod.full, VL_lec_0.mod.red, VL_lec_0.mod.null,  test = "F")
AICctab(VL_lec_0.mod.full, VL_lec_0.mod.red, VL_lec_0.mod.null)


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
lec_1.t1 <- wilcox.test(T1$SES_lec_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_1.t1
lec_1.t8 <- wilcox.test(T8$SES_lec_1, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_1.t8

## compare among treatments and landscape variables
dotchart(SES$SES_lec_1, group = SES$trmt, pch = 19)
with(FS_comps, bartlett.test(SES_lec_1 ~ trmt))
with(KT_comps, bartlett.test(SES_lec_1 ~ trmt))
with(VL_comps, bartlett.test(SES_lec_1 ~ trmt))
with(FS_comps, ad.test(SES_lec_1))
#So normality test is failed. Not sure what to do about that.
#I think that outlier is the cause. It is more than 2 standard devs from the rest of the data points
with(KT_comps, ad.test(SES_lec_1))
with(VL_comps, ad.test(SES_lec_1))

#Okay time to attempt fixing the normality issue (12.27.22)
FS_comps$SES_lec_1 #Did that to see the row with the outlier. Looks like its #1
FSwOUT.out<-FS_comps[-1,]
with(FSwOUT.out, bartlett.test(SES_lec_1 ~ trmt))
with(FSwOUT.out, ad.test(SES_lec_1))
#Okay normality test not failed now. 
#Outlier removed- Time to make a model
FS_lec_1.mod.full<-glm(SES_lec_1 ~ trmt + pland + enn, family = gaussian, data = FSwOUT.out)
summary(FS_lec_1.mod.full)
step(FS_lec_1.mod.full)
#null model is best comparison

FS_lec_1.mod.null <- glm(SES_lec_1 ~ 1, family = gaussian, data = FSwOUT.out)
summary(FS_lec_1.mod.null)
qqnorm(resid(FS_lec_1.mod.null))
qqline(resid(FS_lec_1.mod.null))
plot(simulateResiduals(FS_lec_1.mod.null))
densityPlot(rstudent(FS_lec_1.mod.null)) # check density estimate of the distribution of residuals
outlierTest(FS_lec_1.mod.null)
influenceIndexPlot(FS_lec_1.mod.null, vars = c("Cook"), id = list(n = 3))
#all distances lower than .5 so looks good

# model comparison techniques
anova(FS_lec_1.mod.full, FS_lec_1.mod.null,  test = "F")
AICctab(FS_lec_1.mod.full, FS_lec_1.mod.null)


KT_lec_1.mod.full <- glm(SES_lec_1 ~ trmt + pland + enn, family = gaussian, data = KT_comps)
summary(KT_lec_1.mod.full)
step(KT_lec_1.mod.full)

KT_lec_1.mod.null <- glm(SES_lec_1 ~ 1, family = gaussian, data = KT_comps)
summary(KT_lec_1.mod.null)
qqnorm(resid(KT_lec_1.mod.null))
qqline(resid(KT_lec_1.mod.null))
plot(simulateResiduals(KT_lec_1.mod.null))
densityPlot(rstudent(KT_lec_1.mod.null)) # check density estimate of the distribution of residuals
outlierTest(KT_lec_1.mod.null)
influenceIndexPlot(KT_lec_1.mod.null, vars = c("Cook"), id = list(n = 3))


# model comparison techniques
anova(KT_lec_1.mod.full, KT_lec_1.mod.null,  test = "F")
AICctab(KT_lec_1.mod.full, KT_lec_1.mod.null)


VL_lec_1.mod.full <- glm(SES_lec_1 ~ trmt + pland  + enn, family = gaussian, data = VL_comps)
summary(VL_lec_1.mod.full)
step(VL_lec_1.mod.full)

VL_lec_1.mod.null <- glm(SES_lec_1 ~ 1, family = gaussian, data = VL_comps)
summary(VL_lec_1.mod.null)
qqnorm(resid(VL_lec_1.mod.null))
qqline(resid(VL_lec_1.mod.null))
plot(simulateResiduals(VL_lec_1.mod.null))
densityPlot(rstudent(VL_lec_1.mod.null)) # check density estimate of the distribution of residuals
outlierTest(VL_lec_1.mod.null)
influenceIndexPlot(VL_lec_1.mod.null, vars = c("Cook"), id = list(n = 3))

VL_lec_1.mod.null <- glm(SES_lec_1 ~ 1, family = gaussian, data = VL_comps)

# model comparison techniques
anova(VL_lec_1.mod.full, VL_lec_1.mod.null, VL_lec_1.mod.null,  test = "F")
AICctab(VL_lec_1.mod.full, VL_lec_1.mod.null, VL_lec_1.mod.null)


#I am leaving in the following code because it has outlier removal
#It is useful code to refer back to, just in case
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
lec_2.t1 <- wilcox.test(T1$SES_lec_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_2.t1
lec_2.t8 <- wilcox.test(T8$SES_lec_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
lec_2.t8


## compare among treatments and landscape variables
dotchart(SES$SES_lec_2, group = SES$trmt, pch = 19)
with(FS_comps, bartlett.test(SES_lec_2 ~ trmt))
with(KT_comps, bartlett.test(SES_lec_2 ~ trmt))
with(VL_comps, bartlett.test(SES_lec_2 ~ trmt))
with(FS_comps, ad.test(SES_lec_2))
with(KT_comps, ad.test(SES_lec_2))
with(VL_comps, ad.test(SES_lec_2))


FS_lec_2.mod.full <- glm(SES_lec_2 ~ trmt + pland + enn, family = gaussian, data = FS_comps)
summary(FS_lec_2.mod.full)
step(FS_lec_2.mod.full)

FS_lec_2.mod.red <- glm(SES_lec_2 ~ pland, family = gaussian, data = FS_comps)
qqnorm(resid(FS_lec_2.mod.red))
qqline(resid(FS_lec_2.mod.red))
plot(simulateResiduals(FS_lec_2.mod.red))
densityPlot(rstudent(FS_lec_2.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_lec_2.mod.red)
influenceIndexPlot(FS_lec_2.mod.red, vars = c("Cook"), id = list(n = 3))
#Failed outlier test, but influence is only 0.6. I might just leave it

Anova(FS_lec_2.mod.red)
plot_summs(FS_lec_2.mod.red, scale = TRUE)

effect_plot(FS_lec_2.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = ' (PLAND)', y.label = 'Standardized Effect Sizes (SES)')


FS_lec_2.mod.null <- glm(SES_lec_2 ~ 1, family = gaussian, data = FS_comps)

# model comparison techniques
anova(FS_lec_2.mod.full, FS_lec_2.mod.red, FS_lec_2.mod.null, test = "F")
AICctab(FS_lec_2.mod.full, FS_lec_2.mod.red, FS_lec_2.mod.null)

KT_lec_2.mod.full <- glm(SES_lec_2 ~ trmt + pland + enn, family = gaussian, data = KT_comps)
summary(KT_lec_2.mod.full)
step(KT_lec_2.mod.full)

KT_lec_2.mod.red <- glm(SES_lec_2 ~ pland, family = gaussian, data = KT_comps)
summary(KT_lec_2.mod.red)
qqnorm(resid(KT_lec_2.mod.red))
qqline(resid(KT_lec_2.mod.red))
plot(simulateResiduals(KT_lec_2.mod.red))
densityPlot(rstudent(KT_lec_2.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_lec_2.mod.red)
influenceIndexPlot(KT_lec_2.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(KT_lec_2.mod.red)

plot_summs(KT_lec_2.mod.red, scale = TRUE)

effect_plot(KT_lec_2.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = 'Largest Patch Index (LPI)', y.label = 'Standardized Effect Sizes (SES)')

KT_lec_2.mod.null <- glm(SES_lec_2 ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_lec_2.mod.full, KT_lec_2.mod.red, KT_lec_2.mod.null, test = "F")
AICctab(KT_lec_2.mod.full, KT_lec_2.mod.red, KT_lec_2.mod.null)


VL_lec_2.mod.full <- glm(SES_lec_2 ~ trmt + pland + enn, family = gaussian, data = VL_comps)
summary(VL_lec_2.mod.full)
step(VL_lec_2.mod.full)

VL_lec_2.mod.red <- glm(SES_lec_2 ~ trmt + pland, family = gaussian, data = VL_comps)
summary(VL_lec_2.mod.red)
qqnorm(resid(VL_lec_2.mod.red))
qqline(resid(VL_lec_2.mod.red))
plot(simulateResiduals(VL_lec_2.mod.red))
densityPlot(rstudent(VL_lec_2.mod.red)) # check density estimate of the distribution of residuals
outlierTest(VL_lec_2.mod.red)
influenceIndexPlot(VL_lec_2.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(VL_lec_2.mod.red)
emmeans(VL_lec_2.mod.red, pairwise ~ trmt)
plot_summs(VL_lec_2.mod.red, scale = TRUE)

effect_plot(VL_lec_2.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'TRT', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(VL_lec_2.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = ' (PLAND)', y.label = 'Standardized Effect Sizes (SES)')

VL_lec_2.mod.null <- glm(SES_lec_2 ~ 1, family = gaussian, data = VL_comps)

# model comparison techniques
anova(VL_lec_2.mod.full, VL_lec_2.mod.red, VL_lec_2.mod.null, test = "F")
AICctab(VL_lec_2.mod.full, VL_lec_2.mod.red, VL_lec_2.mod.null)






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
with(FS_comps, bartlett.test(SES_nest_1 ~ trmt))
with(KT_comps, bartlett.test(SES_nest_1 ~ trmt))
with(VL_comps, bartlett.test(SES_nest_1 ~ trmt))
with(FS_comps, ad.test(SES_nest_1))
with(KT_comps, ad.test(SES_nest_1))
with(VL_comps, ad.test(SES_nest_1))


FS_nest_1.mod.full <- glm(SES_nest_1 ~ trmt + pland + enn, family = gaussian, data = FS_comps)
summary(FS_nest_1.mod.full)
step(FS_nest_1.mod.full)

FS_nest_1.mod.red <- glm(SES_nest_1 ~ trmt, family = gaussian, data = FS_comps)
summary(FS_nest_1.mod.red)
qqnorm(resid(FS_nest_1.mod.red))
qqline(resid(FS_nest_1.mod.red))
plot(simulateResiduals(FS_nest_1.mod.red))
densityPlot(rstudent(FS_nest_1.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_nest_1.mod.red)
influenceIndexPlot(FS_nest_1.mod.red, vars = c("Cook"), id = list(n = 3))

#it looks like there is an outlier
FS_nest_1.mod.red2 <- update(FS_nest_1.mod.red, subset = -c(1))
summary(FS_nest_1.mod.red2)
compareCoefs(FS_nest_1.mod.red, FS_nest_1.mod.red2) # compares estimated coefficients and their standard errors

Anova(FS_nest_1.mod.red)
emmeans(FS_nest_1.mod.red, pairwise ~ trmt)

effect_plot(FS_nest_1.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')

FS_nest_1.mod.null <- glm(SES_nest_1 ~ 1, family = gaussian, data = FS_comps)

# model comparison techniques
anova(FS_nest_1.mod.full, FS_nest_1.mod.red, FS_nest_1.mod.null, test = "F")
AICctab(FS_nest_1.mod.full, FS_nest_1.mod.red, FS_nest_1.mod.null)


KT_nest_1.mod.full <- glm(SES_nest_1 ~ trmt + pland + enn, family = gaussian, data = KT_comps)
summary(KT_nest_1.mod.full)
step(KT_nest_1.mod.full)

KT_nest_1.mod.red <- glm(SES_nest_1 ~ trmt+pland, family = gaussian, data =  KT_comps)
summary(KT_nest_1.mod.red)
qqnorm(resid(KT_nest_1.mod.red))
qqline(resid(KT_nest_1.mod.red))
plot(simulateResiduals(KT_nest_1.mod.red))
densityPlot(rstudent(KT_nest_1.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_nest_1.mod.red)
influenceIndexPlot(KT_nest_1.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(KT_nest_1.mod.red)
emmeans(KT_nest_1.mod.red, pairwise ~ trmt)

effect_plot(KT_nest_1.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(KT_nest_1.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')


KT_nest_1.mod.null <- glm(SES_nest_1 ~ 1, family = gaussian, data = 
                            KT_comps)

# model comparison techniques
anova(KT_nest_1.mod.full, KT_nest_1.mod.red, KT_nest_1.mod.null, test = "F")
AICctab(KT_nest_1.mod.full, KT_nest_1.mod.red, KT_nest_1.mod.null)

#Below is kruskal test for nonparametric data (Thank you Kayla)
kruskal.test(SES_nest_1 ~ trmt, data = VL_comps)
#Basically for these violations of assumption, because it is between our different study vacant lots, we are just looking at significance of trt



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
nest_2.t1 <- wilcox.test(T1$SES_nest_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_2.t1
nest_2.t8 <- wilcox.test(T8$SES_nest_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_2.t8


## compare among treatments and landscape variables
dotchart(SES$SES_nest_2, group = SES$trmt, pch = 19)
with(FS_comps, bartlett.test(SES_nest_2 ~ trmt))
with(KT_comps, bartlett.test(SES_nest_2 ~ trmt))
with(VL_comps, bartlett.test(SES_nest_2 ~ trmt))
with(FS_comps, ad.test(SES_nest_2))
with(KT_comps, ad.test(SES_nest_2))
with(VL_comps, ad.test(SES_nest_2))


FS_nest_2.mod.full <- glm(SES_nest_2 ~ trmt + pland + enn, family = gaussian, data = FS_comps)
summary(FS_nest_2.mod.full)
step(FS_nest_2.mod.full)

FS_nest_2.mod.red <- glm(SES_nest_2 ~ trmt, family = gaussian, data = FS_comps)
summary(FS_nest_2.mod.red)
qqnorm(resid(FS_nest_2.mod.red))
qqline(resid(FS_nest_2.mod.red))
plot(simulateResiduals(FS_nest_2.mod.red))
densityPlot(rstudent(FS_nest_2.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_nest_2.mod.red)
influenceIndexPlot(FS_nest_2.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(FS_nest_2.mod.red)
emmeans(FS_nest_2.mod.red, pairwise ~ trmt)

effect_plot(FS_nest_2.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')

FS_nest_2.mod.null <- glm(SES_nest_2 ~ 1, family = gaussian, data = FS_comps)

# model comparison techniques
anova(FS_nest_2.mod.full, FS_nest_2.mod.red, FS_nest_2.mod.null, test = "F")
AICctab(FS_nest_2.mod.full, FS_nest_2.mod.red, FS_nest_2.mod.null)


KT_nest_2.mod.full <- glm(SES_nest_2 ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_nest_2.mod.full)
step(KT_nest_2.mod.full)

KT_nest_2.mod.null <- glm(SES_nest_2 ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_nest_2.mod.full, KT_nest_2.mod.null, test = "F")
AICctab(KT_nest_2.mod.full, KT_nest_2.mod.null)


VL_nest_2.mod.full <- glm(SES_nest_2 ~ trmt + pland + enn, family = gaussian, data = VL_comps)
summary(VL_nest_2.mod.full)
step(VL_nest_2.mod.full)


VL_nest_2.mod.red <- glm(SES_nest_2 ~ trmt, family = gaussian, data = VL_comps)
summary(VL_nest_2.mod.red)
qqnorm(resid(VL_nest_2.mod.red))
qqline(resid(VL_nest_2.mod.red))
plot(simulateResiduals(VL_nest_2.mod.red))
densityPlot(rstudent(VL_nest_2.mod.red)) # check density estimate of the distribution of residuals
outlierTest(VL_nest_2.mod.red)
influenceIndexPlot(VL_nest_2.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(VL_nest_2.mod.red)
emmeans(VL_nest_2.mod.red, pairwise ~ trmt)

effect_plot(VL_nest_2.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')

VL_nest_2.mod.null <- glm(SES_nest_2 ~ 1, family = gaussian, data = VL_comps)
# model comparison techniques
anova(VL_nest_2.mod.red, VL_nest_2.mod.full, VL_nest_2.mod.null, test = "F")
AICctab(VL_nest_2.mod.red, VL_nest_2.mod.full, VL_nest_2.mod.null)



### nest_3 - Hive----
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
with(FS_comps, bartlett.test(SES_nest_3 ~ trmt))
with(KT_comps, bartlett.test(SES_nest_3 ~ trmt))
with(VL_comps, bartlett.test(SES_nest_3 ~ trmt))
with(FS_comps, ad.test(SES_nest_3))
with(KT_comps, ad.test(SES_nest_3))
with(VL_comps, ad.test(SES_nest_3))


FS_nest_3.mod.full <- glm(SES_nest_3 ~ trmt + pland  + enn, family = gaussian, data = FS_comps)
summary(FS_nest_3.mod.full)
step(FS_nest_3.mod.full)

FS_nest_3.mod.red <- glm(SES_nest_3 ~ pland+enn, family = gaussian, data = FS_comps)
summary(FS_nest_3.mod.red)
qqnorm(resid(FS_nest_3.mod.red))
qqline(resid(FS_nest_3.mod.red))
plot(simulateResiduals(FS_nest_3.mod.red))
#I am adding the below code to look at the graphs individually. the second one has red lines
plotQQunif(FS_nest_3.mod.red)
plotResiduals(FS_nest_3.mod.red)
densityPlot(rstudent(FS_nest_3.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_nest_3.mod.red)
influenceIndexPlot(FS_nest_3.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(FS_nest_3.mod.red)

effect_plot(FS_nest_3.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(FS_nest_3.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')


FS_nest_3.mod.null <- glm(SES_nest_3 ~ 1, family = gaussian, data = FS_comps)

# model comparison techniques
anova(FS_nest_3.mod.full, FS_nest_3.mod.red, FS_nest_3.mod.null, test = "F")
AICctab(FS_nest_3.mod.full, FS_nest_3.mod.red, FS_nest_3.mod.null)


KT_nest_3.mod.full <- glm(SES_nest_3 ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_nest_3.mod.full)
step(KT_nest_3.mod.full)

KT_nest_3.mod.red <- glm(SES_nest_3 ~ trmt, family = gaussian, data = KT_comps)
summary(KT_nest_3.mod.red)
qqnorm(resid(KT_nest_3.mod.red))
qqline(resid(KT_nest_3.mod.red))
plot(simulateResiduals(KT_nest_3.mod.red))
densityPlot(rstudent(KT_nest_3.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_nest_3.mod.red)
influenceIndexPlot(KT_nest_3.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(KT_nest_3.mod.red)
emmeans(KT_nest_3.mod.red, pairwise ~ trmt)

effect_plot(KT_nest_3.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')


KT_nest_3.mod.null <- glm(SES_nest_3 ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_nest_3.mod.full, KT_nest_3.mod.red, KT_nest_3.mod.null, test = "F")
AICctab(KT_nest_3.mod.full, KT_nest_3.mod.red, KT_nest_3.mod.null)


VL_nest_3.mod.full <- glm(SES_nest_3 ~ trmt + pland  + enn, family = gaussian, data = VL_comps)
summary(VL_nest_3.mod.full)
step(VL_nest_3.mod.full)

VL_nest_3.mod.null <- glm(SES_nest_3 ~ 1, family = gaussian, data = VL_comps)
summary(VL_nest_3.mod.null)
qqnorm(resid(VL_nest_3.mod.null))
qqline(resid(VL_nest_3.mod.null))
plot(simulateResiduals(VL_nest_3.mod.null))
densityPlot(rstudent(VL_nest_3.mod.null)) # check density estimate of the distribution of residuals
outlierTest(VL_nest_3.mod.null)
influenceIndexPlot(VL_nest_3.mod.null, vars = c("Cook"), id = list(n = 3))


# model comparison techniques
anova(VL_nest_3.mod.full, VL_nest_3.mod.null, test = "F")
AICctab(VL_nest_3.mod.full, VL_nest_3.mod.null)






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
nest_4.t1 <- wilcox.test(T1$SES_nest_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_4.t1
nest_4.t8 <- wilcox.test(T8$SES_nest_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_4.t8


## compare among treatments and landscape variables
dotchart(SES$SES_nest_4, group = SES$trmt, pch = 19)
with(FS_comps, bartlett.test(SES_nest_4 ~ trmt))
with(KT_comps, bartlett.test(SES_nest_4 ~ trmt))
with(VL_comps, bartlett.test(SES_nest_4 ~ trmt))
#Violation for VL- KW test will be used instead
with(FS_comps, ad.test(SES_nest_4))
with(KT_comps, ad.test(SES_nest_4))
with(VL_comps, ad.test(SES_nest_4))


FS_nest_4.mod.full <- glm(SES_nest_4 ~ trmt + pland  + enn, family = gaussian, data = FS_comps)
summary(FS_nest_4.mod.full)
step(FS_nest_4.mod.full)

FS_nest_4.mod.null <- glm(SES_nest_4 ~ 1, family = gaussian, data = FS_comps)
summary(FS_nest_4.mod.null)
qqnorm(resid(FS_nest_4.mod.null))
qqline(resid(FS_nest_4.mod.null))
plot(simulateResiduals(FS_nest_4.mod.null))
densityPlot(rstudent(FS_nest_4.mod.null)) # check density estimate of the distribution of residuals
outlierTest(FS_nest_4.mod.null)
influenceIndexPlot(FS_nest_4.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(FS_nest_4.mod.full,  FS_nest_4.mod.null, test = "F")
AICctab(FS_nest_4.mod.full, FS_nest_4.mod.null)


KT_nest_4.mod.full <- glm(SES_nest_4 ~ trmt + pland + enn, family = gaussian, data = KT_comps)
summary(KT_nest_4.mod.full)
step(KT_nest_4.mod.full)

KT_nest_4.mod.red <- glm(SES_nest_4 ~ trmt+pland+enn, family = gaussian, data = KT_comps)
summary(KT_nest_4.mod.red)
qqnorm(resid(KT_nest_4.mod.red))
qqline(resid(KT_nest_4.mod.red))
plot(simulateResiduals(KT_nest_4.mod.red))
densityPlot(rstudent(KT_nest_4.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_nest_4.mod.red)
influenceIndexPlot(KT_nest_4.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(KT_nest_4.mod.red)
emmeans(KT_nest_4.mod.red, pairwise ~ trmt)

effect_plot(KT_nest_4.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'TRT', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(KT_nest_4.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = 'PLAND', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(KT_nest_4.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'ENN', y.label = 'Standardized Effect Sizes (SES)')


KT_nest_4.mod.null <- glm(SES_nest_4 ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_nest_4.mod.full, KT_nest_4.mod.red, KT_nest_4.mod.null, test = "F")
AICctab(KT_nest_4.mod.full, KT_nest_4.mod.red, KT_nest_4.mod.null)
#The full model and reduced model are literally the same.

kruskal.test(SES_nest_4 ~ trmt, data = VL_comps)#for the vacant lot comparison, you can use a kruskal-wallis test
#if the data violate the normality or homogeneity assumptions. Kruskal-wallis is a non-parametric test
#thank you to Kayla for the above code

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
nest_5.t1 <- wilcox.test(T1$SES_nest_5, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_5.t1
nest_5.t8 <- wilcox.test(T8$SES_nest_5, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
nest_5.t8


## compare among treatments and landscape variables
dotchart(SES$SES_nest_5, group = SES$trmt, pch = 19)
with(FS_comps, bartlett.test(SES_nest_5 ~ trmt))
with(KT_comps, bartlett.test(SES_nest_5 ~ trmt))
with(VL_comps, bartlett.test(SES_nest_5 ~ trmt))
with(FS_comps, ad.test(SES_nest_5))
with(KT_comps, ad.test(SES_nest_5))
with(VL_comps, ad.test(SES_nest_5))


FS_nest_5.mod.full <- glm(SES_nest_5 ~ trmt + pland  + enn, family = gaussian, data = FS_comps)
summary(FS_nest_5.mod.full)
step(FS_nest_5.mod.full)

FS_nest_5.mod.null <- glm(SES_nest_5 ~ 1, family = gaussian, data = FS_comps)
summary(FS_nest_5.mod.null)
qqnorm(resid(FS_nest_5.mod.null))
qqline(resid(FS_nest_5.mod.null))
plot(simulateResiduals(FS_nest_5.mod.null))
densityPlot(rstudent(FS_nest_5.mod.null)) # check density estimate of the distribution of residuals
outlierTest(FS_nest_5.mod.null)
influenceIndexPlot(FS_nest_5.mod.null, vars = c("Cook"), id = list(n = 3))


# model comparison techniques
anova(FS_nest_5.mod.full,  FS_nest_5.mod.null, test = "F")
AICctab(FS_nest_5.mod.full,  FS_nest_5.mod.null)


KT_nest_5.mod.full <- glm(SES_nest_5 ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_nest_5.mod.full)
step(KT_nest_5.mod.full)

KT_nest_5.mod.red <- glm(SES_nest_5 ~ trmt, family = gaussian, data = KT_comps)
summary(KT_nest_5.mod.red)
qqnorm(resid(KT_nest_5.mod.red))
qqline(resid(KT_nest_5.mod.red))
plot(simulateResiduals(KT_nest_5.mod.red))
densityPlot(rstudent(KT_nest_5.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_nest_5.mod.red)
influenceIndexPlot(KT_nest_5.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(KT_nest_5.mod.red)
emmeans(KT_nest_5.mod.red, pairwise ~ trmt)

effect_plot(KT_nest_5.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')
#for some reason I keep needing to do dev.off() in order to get the effect plots to run. but it works so yaY?

KT_nest_5.mod.null <- glm(SES_nest_5 ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_nest_5.mod.full, KT_nest_5.mod.red, KT_nest_5.mod.null, test = "F")
AICctab(KT_nest_5.mod.full, KT_nest_5.mod.red, KT_nest_5.mod.null)

VL_nest_5.mod.full <- glm(SES_nest_5 ~ trmt + pland + enn, family = gaussian, data = VL_comps)
summary(VL_nest_5.mod.full)
step(VL_nest_5.mod.full)

VL_nest_5.mod.red <- glm(SES_nest_5 ~ enn, family = gaussian, data = VL_comps)
summary(VL_nest_5.mod.red)
qqnorm(resid(VL_nest_5.mod.red))
qqline(resid(VL_nest_5.mod.red))
plot(simulateResiduals(VL_nest_5.mod.red))
densityPlot(rstudent(VL_nest_5.mod.red)) # check density estimate of the distribution of residuals
outlierTest(VL_nest_5.mod.red)
influenceIndexPlot(VL_nest_5.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(VL_nest_5.mod.red)

effect_plot(VL_nest_5.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'ENN', y.label = 'Standardized Effect Sizes (SES)')


VL_nest_5.mod.null <- glm(SES_nest_5 ~ 1, family = gaussian, data = VL_comps)

# model comparison techniques
anova(VL_nest_5.mod.full, VL_nest_5.mod.red, VL_nest_5.mod.null, test = "F")
AICctab(VL_nest_5.mod.full, VL_nest_5.mod.red, VL_nest_5.mod.null)



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

with(FS_comps, bartlett.test(SES_soc_1 ~ trmt))
with(KT_comps, bartlett.test(SES_soc_1 ~ trmt))
with(VL_comps, bartlett.test(SES_soc_1 ~ trmt))
with(FS_comps, ad.test(SES_soc_1))
with(KT_comps, ad.test(SES_soc_1))
with(VL_comps, ad.test(SES_soc_1))


FS_soc_1.mod.full <- glm(SES_soc_1 ~ trmt + pland  + enn, family = gaussian, data = FS_comps)
summary(FS_soc_1.mod.full)
step(FS_soc_1.mod.full)

FS_soc_1.mod.red <- glm(SES_soc_1 ~ trmt+enn, family = gaussian, data = FS_comps)
summary(FS_soc_1.mod.red)
qqnorm(resid(FS_soc_1.mod.red))
qqline(resid(FS_soc_1.mod.red))
plot(simulateResiduals(FS_soc_1.mod.red))
densityPlot(rstudent(FS_soc_1.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_soc_1.mod.red)
influenceIndexPlot(FS_soc_1.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(FS_soc_1.mod.red)
emmeans(FS_soc_1.mod.red, pairwise ~ trmt)

dev.off()
effect_plot(FS_soc_1.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(FS_soc_1.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')


FS_soc_1.mod.null <- glm(SES_soc_1 ~ 1, family = gaussian, data = FS_comps)

# model comparison techniques
anova(FS_soc_1.mod.full, FS_soc_1.mod.red, FS_soc_1.mod.null, test = "F")
AICctab(FS_soc_1.mod.full, FS_soc_1.mod.red, FS_soc_1.mod.null)


KT_soc_1.mod.full <- glm(SES_soc_1 ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_soc_1.mod.full)
step(KT_soc_1.mod.full)

KT_soc_1.mod.red <- glm(SES_soc_1 ~ pland+ enn, family = gaussian, data = KT_comps)
summary(KT_soc_1.mod.red)
qqnorm(resid(KT_soc_1.mod.red))
qqline(resid(KT_soc_1.mod.red))
plot(simulateResiduals(KT_soc_1.mod.red))
densityPlot(rstudent(KT_soc_1.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_soc_1.mod.red)
influenceIndexPlot(KT_soc_1.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(KT_soc_1.mod.red)

dev.off()
effect_plot(KT_soc_1.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = 'PLAND', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(KT_soc_1.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'ENN', y.label = 'Standardized Effect Sizes (SES)')


KT_soc_1.mod.null <- glm(SES_soc_1 ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_soc_1.mod.full, KT_soc_1.mod.red, KT_soc_1.mod.null, test = "F")
AICctab(KT_soc_1.mod.full, KT_soc_1.mod.red, KT_soc_1.mod.null)

VL_soc_1.mod.full <- glm(SES_soc_1 ~ trmt + pland + enn, family = gaussian, data = VL_comps)
summary(VL_soc_1.mod.full)
step(VL_soc_1.mod.full)

VL_soc_1.mod.red <- glm(SES_soc_1 ~ trmt+pland+enn, family = gaussian, data = VL_comps)
summary(VL_soc_1.mod.red)
qqnorm(resid(VL_soc_1.mod.red))
qqline(resid(VL_soc_1.mod.red))
plot(simulateResiduals(VL_soc_1.mod.red))
densityPlot(rstudent(VL_soc_1.mod.red)) # check density estimate of the distribution of residuals
outlierTest(VL_soc_1.mod.red)
influenceIndexPlot(VL_soc_1.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(VL_soc_1.mod.red)
emmeans(VL_soc_1.mod.red, pairwise ~ trmt)

dev.off()
effect_plot(VL_soc_1.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'trt', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(VL_soc_1.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = 'pland', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(VL_soc_1.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'enn', y.label = 'Standardized Effect Sizes (SES)')


VL_soc_1.mod.null <- glm(SES_soc_1 ~ 1, family = gaussian, data = VL_comps)

# model comparison techniques
anova(VL_soc_1.mod.full, VL_soc_1.mod.red, VL_soc_1.mod.null, test = "F")
AICctab(VL_soc_1.mod.full, VL_soc_1.mod.red, VL_soc_1.mod.null)



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
soc_2.t1 <- wilcox.test(T1$SES_soc_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_2.t1
soc_2.t8 <- wilcox.test(T8$SES_soc_2, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_2.t8


## compare among treatments and landscape variables
dotchart(SES$SES_soc_2, group = SES$trmt, pch = 19)

with(FS_comps, bartlett.test(SES_soc_2 ~ trmt))
with(KT_comps, bartlett.test(SES_soc_2 ~ trmt))
with(VL_comps, bartlett.test(SES_soc_2 ~ trmt))
with(FS_comps, ad.test(SES_soc_2))
with(KT_comps, ad.test(SES_soc_2))
with(VL_comps, ad.test(SES_soc_2))
#Normality violated- Are there any outliers though?
VL_comps$SES_soc_2 
mean(VL_comps$SES_soc_2)
sd(VL_comps$SES_soc_2)
-0.8309074+(3*1.116523)
#nope. looks like data point 9 is close, but not technically an outlier
#so we will do KW test

FS_soc_2.mod.full <- glm(SES_soc_2 ~ trmt + pland  + enn, family = gaussian, data = FS_comps)
summary(FS_soc_2.mod.full)
step(FS_soc_2.mod.full)

FS_soc_2.mod.null <- glm(SES_soc_2 ~ 1, family = gaussian, data = FS_comps)
summary(FS_soc_2.mod.null)
qqnorm(resid(FS_soc_2.mod.null))
qqline(resid(FS_soc_2.mod.null))
plot(simulateResiduals(FS_soc_2.mod.null))
densityPlot(rstudent(FS_soc_2.mod.null)) # check density estimate of the distribution of residuals
outlierTest(FS_soc_2.mod.null)
influenceIndexPlot(FS_soc_2.mod.null, vars = c("Cook"), id = list(n = 3))
#significant outlier, let's remove it and see if the model fits better
FS_soc_2.mod.red <- update(FS_soc_2.mod.null, subset = -c(13))
summary(FS_soc_2.mod.red)
compareCoefs(FS_soc_2.mod.null, FS_soc_2.mod.red)

#out of curiousity- is the outlier in the original model?
outlierTest(FS_soc_2.mod.full)
influenceIndexPlot(FS_soc_2.mod.full, vars = c("Cook"), id = list(n = 3))
#yep- okay just trying the following to sate my own curiosity
FS_soc_2.mod.full2 <- update(FS_soc_2.mod.full, subset = -c(11))
FS_soc_2.mod.full2 <- update(FS_soc_2.mod.full2, subset = -c(13))
summary(FS_soc_2.mod.full2)
step(FS_soc_2.mod.full2)
#okay, that makes intercept the best fit


# model comparison techniques
anova(FS_soc_2.mod.full,  FS_soc_2.mod.null, test = "F")
AICctab(FS_soc_2.mod.full,  FS_soc_2.mod.null)



KT_soc_2.mod.full <- glm(SES_soc_2 ~ trmt + pland + enn, family = gaussian, data = KT_comps)
summary(KT_soc_2.mod.full)
step(KT_soc_2.mod.full)

KT_soc_2.mod.null <- glm(SES_soc_2 ~ 1, family = gaussian, data = KT_comps)
summary(KT_soc_2.mod.null)
qqnorm(resid(KT_soc_2.mod.null))
qqline(resid(KT_soc_2.mod.null))
plot(simulateResiduals(KT_soc_2.mod.null))
densityPlot(rstudent(KT_soc_2.mod.null)) # check density estimate of the distribution of residuals
outlierTest(KT_soc_2.mod.null)
influenceIndexPlot(KT_soc_2.mod.null, vars = c("Cook"), id = list(n = 3))


# model comparison techniques
anova(KT_soc_2.mod.full, KT_soc_2.mod.null)
AICctab(KT_soc_2.mod.full, KT_soc_2.mod.null)


#Just testing trt because normality violated
kruskal.test(SES_soc_2 ~ trmt, data = VL_comps)





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
soc_3.t1 <- wilcox.test(T1$SES_soc_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_3.t1
soc_3.t8 <- wilcox.test(T8$SES_soc_3, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_3.t8


## compare among treatments and landscape variables
dotchart(SES$SES_soc_3, group = SES$trmt, pch = 19)
with(FS_comps, bartlett.test(SES_soc_3 ~ trmt))
with(KT_comps, bartlett.test(SES_soc_3 ~ trmt))
with(VL_comps, bartlett.test(SES_soc_3 ~ trmt))
with(FS_comps, ad.test(SES_soc_3))
with(KT_comps, ad.test(SES_soc_3))
with(VL_comps, ad.test(SES_soc_3))


FS_soc_3.mod.full <- glm(SES_soc_3 ~ trmt + pland  + enn, family = gaussian, data = FS_comps)
summary(FS_soc_3.mod.full)
step(FS_soc_3.mod.full)

FS_soc_3.mod.red <- glm(SES_soc_3 ~ trmt+pland+enn, family = gaussian, data = FS_comps)
summary(FS_soc_3.mod.red)
qqnorm(resid(FS_soc_3.mod.red))
qqline(resid(FS_soc_3.mod.red))
plot(simulateResiduals(FS_soc_3.mod.red))
densityPlot(rstudent(FS_soc_3.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_soc_3.mod.red)
influenceIndexPlot(FS_soc_3.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(FS_soc_3.mod.red)

effect_plot(FS_soc_3.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(FS_soc_3.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')
effect_plot(FS_soc_3.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')


FS_soc_3.mod.null <- glm(SES_soc_3 ~ 1, family = gaussian, data = FS_comps)

# model comparison techniques
anova(FS_soc_3.mod.full, FS_soc_3.mod.red, FS_soc_3.mod.null, test = "F")
AICctab(FS_soc_3.mod.full, FS_soc_3.mod.red, FS_soc_3.mod.null)


KT_soc_3.mod.full <- glm(SES_soc_3 ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_soc_3.mod.full)
step(KT_soc_3.mod.full)

KT_soc_3.mod.red <- glm(SES_soc_3 ~ pland, family = gaussian, data = KT_comps)
summary(KT_soc_3.mod.red)
qqnorm(resid(KT_soc_3.mod.red))
qqline(resid(KT_soc_3.mod.red))
plot(simulateResiduals(KT_soc_3.mod.red))
densityPlot(rstudent(KT_soc_3.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_soc_3.mod.red)
influenceIndexPlot(KT_soc_3.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(KT_soc_3.mod.red)

effect_plot(KT_soc_3.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = 'Percentage Greenspace (PLAND)', y.label = 'Standardized Effect Sizes (SES)')


KT_soc_3.mod.null <- glm(SES_soc_3 ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_soc_3.mod.full, KT_soc_3.mod.red, KT_soc_3.mod.null, test = "F")
AICctab(KT_soc_3.mod.full, KT_soc_3.mod.red, KT_soc_3.mod.null)


VL_soc_3.mod.full <- glm(SES_soc_3 ~ trmt + pland  + enn, family = gaussian, data = VL_comps)
summary(VL_soc_3.mod.full)
step(VL_soc_3.mod.full)

VL_soc_3.mod.null <- glm(SES_soc_3 ~ 1, family = gaussian, data = VL_comps)
summary(VL_soc_3.mod.null)
qqnorm(resid(VL_soc_3.mod.null))
qqline(resid(VL_soc_3.mod.null))
plot(simulateResiduals(VL_soc_3.mod.null))
densityPlot(rstudent(VL_soc_3.mod.null)) # check density estimate of the distribution of residuals
outlierTest(VL_soc_3.mod.null)
influenceIndexPlot(VL_soc_3.mod.null, vars = c("Cook"), id = list(n = 3))


# model comparison techniques
anova(VL_soc_3.mod.full, VL_soc_3.mod.null, test = "F")
AICctab(VL_soc_3.mod.full, VL_soc_3.mod.null)




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
soc_4.t1 <- wilcox.test(T1$SES_soc_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_4.t1
soc_4.t8 <- wilcox.test(T8$SES_soc_4, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
soc_4.t8


## compare among treatments and landscape variables
dotchart(SES$SES_soc_4, group = SES$trmt, pch = 19)
with(FS_comps, bartlett.test(SES_soc_4 ~ trmt))
with(KT_comps, bartlett.test(SES_soc_4 ~ trmt))
with(VL_comps, bartlett.test(SES_soc_4 ~ trmt))
with(FS_comps, ad.test(SES_soc_4))
with(KT_comps, ad.test(SES_soc_4))
with(VL_comps, ad.test(SES_soc_4))


FS_soc_4.mod.full <- glm(SES_soc_4 ~ trmt + pland  + enn, family = gaussian, data = FS_comps)
summary(FS_soc_4.mod.full)
step(FS_soc_4.mod.full)

FS_soc_4.mod.null <- glm(SES_soc_4 ~ 1, family = gaussian, data = FS_comps)
summary(FS_soc_4.mod.null)
qqnorm(resid(FS_soc_4.mod.null))
qqline(resid(FS_soc_4.mod.null))
plot(simulateResiduals(FS_soc_4.mod.null))
densityPlot(rstudent(FS_soc_4.mod.null)) # check density estimate of the distribution of residuals
outlierTest(FS_soc_4.mod.null)
influenceIndexPlot(FS_soc_4.mod.null, vars = c("Cook"), id = list(n = 3))


# model comparison techniques
anova(FS_soc_4.mod.full, FS_soc_4.mod.null, test = "F")
AICctab(FS_soc_4.mod.full, FS_soc_4.mod.null)


KT_soc_4.mod.full <- glm(SES_soc_4 ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_soc_4.mod.full)
step(KT_soc_4.mod.full)

KT_soc_4.mod.red <- glm(SES_soc_4 ~ enn, family = gaussian, data = KT_comps)
summary(KT_soc_4.mod.red)
qqnorm(resid(KT_soc_4.mod.red))
qqline(resid(KT_soc_4.mod.red))
plot(simulateResiduals(KT_soc_4.mod.red))
densityPlot(rstudent(KT_soc_4.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_soc_4.mod.red)
influenceIndexPlot(KT_soc_4.mod.red, vars = c("Cook"), id = list(n = 3))
#okay there is a outlier
KT_soc_4.mod.red2 <- update(KT_soc_4.mod.red, subset = -c(3))
summary(KT_soc_4.mod.red2)
compareCoefs(KT_soc_4.mod.red2, KT_soc_4.mod.red) # compares estimated coefficients and their standard errors



Anova(KT_soc_4.mod.red)

effect_plot(KT_soc_4.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'ENN', y.label = 'Standardized Effect Sizes (SES)')


KT_soc_4.mod.null <- glm(SES_soc_4 ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_soc_4.mod.full, KT_soc_4.mod.red, KT_soc_4.mod.null, test = "F")
AICctab(KT_soc_4.mod.full, KT_soc_4.mod.red, KT_soc_4.mod.null)


VL_soc_4.mod.full <- glm(SES_soc_4 ~ trmt + pland + enn, family = gaussian, data = VL_comps)
summary(VL_soc_4.mod.full)
step(VL_soc_4.mod.full)

VL_soc_4.mod.null <- glm(SES_soc_4 ~ 1, family = gaussian, data = VL_comps)
summary(VL_soc_4.mod.null)
qqnorm(resid(VL_soc_4.mod.null))
qqline(resid(VL_soc_4.mod.null))
plot(simulateResiduals(VL_soc_4.mod.null))
densityPlot(rstudent(VL_soc_4.mod.null)) # check density estimate of the distribution of residuals
outlierTest(VL_soc_4.mod.null)
influenceIndexPlot(VL_soc_4.mod.null, vars = c("Cook"), id = list(n = 3))



# model comparison techniques
anova(VL_soc_4.mod.full,  VL_soc_4.mod.null, test = "F")
AICctab(VL_soc_4.mod.full, VL_soc_4.mod.null)


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
with(FS_comps, bartlett.test(SES_ori_0 ~ trmt))
with(KT_comps, bartlett.test(SES_ori_0 ~ trmt))
with(VL_comps, bartlett.test(SES_ori_0 ~ trmt))
with(FS_comps, ad.test(SES_ori_0))
with(KT_comps, ad.test(SES_ori_0))
with(VL_comps, ad.test(SES_ori_0))


FS_ori_0.mod.full <- glm(SES_ori_0 ~ trmt + pland + enn, family = gaussian, data = FS_comps)
summary(FS_ori_0.mod.full)
step(FS_ori_0.mod.full)


FS_ori_0.mod.red <- glm(SES_ori_0 ~ enn, family = gaussian, data = FS_comps)
summary(FS_ori_0.mod.red)
qqnorm(resid(FS_ori_0.mod.red))
qqline(resid(FS_ori_0.mod.red))
plot(simulateResiduals(FS_ori_0.mod.red))
densityPlot(rstudent(FS_ori_0.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_ori_0.mod.red)
influenceIndexPlot(FS_ori_0.mod.red, vars = c("Cook"), id = list(n = 3))
#okay there is an outlier
FS_ori_0.mod.red2 <- update(FS_ori_0.mod.red, subset = -c(13))
summary(FS_ori_0.mod.red2)
compareCoefs(FS_ori_0.mod.red2, FS_ori_0.mod.red) # compares estimated coefficients and their standard errors
#It doesn't look like they change much

Anova(FS_ori_0.mod.red)


effect_plot(FS_ori_0.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = ' (ENN)', y.label = 'Standardized Effect Sizes (SES)')

FS_ori_0.mod.null <- glm(SES_ori_0 ~ 1, family = gaussian, data = FS_comps)


# model comparison techniques
anova(FS_ori_0.mod.full, FS_ori_0.mod.red, FS_ori_0.mod.null)
AICctab(FS_ori_0.mod.full, FS_ori_0.mod.red, FS_ori_0.mod.null)



KT_ori_0.mod.full <- glm(SES_ori_0 ~ trmt + pland + enn, family = gaussian, data = KT_comps)
summary(KT_ori_0.mod.full)
step(KT_ori_0.mod.full)

KT_ori_0.mod.null <- glm(SES_ori_0 ~ 1, family = gaussian, data = KT_comps)
summary(KT_ori_0.mod.null)
qqnorm(resid(KT_ori_0.mod.null))
qqline(resid(KT_ori_0.mod.null))
plot(simulateResiduals(KT_ori_0.mod.null))
densityPlot(rstudent(KT_ori_0.mod.null)) # check density estimate of the distribution of residuals
outlierTest(KT_ori_0.mod.null)
influenceIndexPlot(KT_ori_0.mod.null, vars = c("Cook"), id = list(n = 3))


KT_ori_0.mod.null <- glm(SES_ori_0 ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_ori_0.mod.full, KT_ori_0.mod.null)
AICctab(KT_ori_0.mod.full,  KT_ori_0.mod.null)


VL_ori_0.mod.full <- glm(SES_ori_0 ~ trmt + pland  + enn, family = gaussian, data = VL_comps)
summary(VL_ori_0.mod.full)
step(VL_ori_0.mod.full)


VL_ori_0.mod.red <- glm(SES_ori_0 ~ enn, family = gaussian, data = VL_comps)
summary(VL_ori_0.mod.red)
qqnorm(resid(VL_ori_0.mod.red))
qqline(resid(VL_ori_0.mod.red))
plot(simulateResiduals(VL_ori_0.mod.red))
densityPlot(rstudent(VL_ori_0.mod.red)) # check density estimate of the distribution of residuals
outlierTest(VL_ori_0.mod.red)
influenceIndexPlot(VL_ori_0.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(VL_ori_0.mod.red)

effect_plot(VL_ori_0.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Largest Patch Index (LPI)', y.label = 'Standardized Effect Sizes (SES)')

VL_ori_0.mod.null <- glm(SES_ori_0 ~ 1, family = gaussian, data = VL_comps)

# model comparison techniques
anova(VL_ori_0.mod.full, VL_ori_0.mod.red, VL_ori_0.mod.null)
AICctab(VL_ori_0.mod.full, VL_ori_0.mod.red, VL_ori_0.mod.null)




### ori_2 - Exotic----
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
with(FS_comps, bartlett.test(SES_ori_1 ~ trmt))
with(KT_comps, bartlett.test(SES_ori_1 ~ trmt))
with(VL_comps, bartlett.test(SES_ori_1 ~ trmt))
with(FS_comps, ad.test(SES_ori_1))
with(KT_comps, ad.test(SES_ori_1))
with(VL_comps, ad.test(SES_ori_1))


FS_ori_1.mod.full <- glm(SES_ori_1 ~ trmt + pland + enn, family = gaussian, data = FS_comps)
summary(FS_ori_1.mod.full)
step(FS_ori_1.mod.full)


FS_ori_1.mod.red <- glm(SES_ori_1 ~ enn, family = gaussian, data = FS_comps)
summary(FS_ori_1.mod.red)
qqnorm(resid(FS_ori_1.mod.red))
qqline(resid(FS_ori_1.mod.red))
plot(simulateResiduals(FS_ori_1.mod.red))
densityPlot(rstudent(FS_ori_1.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_ori_1.mod.red)
influenceIndexPlot(FS_ori_1.mod.red, vars = c("Cook"), id = list(n = 3))
#okay there is an outlier
FS_ori_1.mod.red2 <- update(FS_ori_1.mod.red, subset = -c(13))
summary(FS_ori_1.mod.red2)
compareCoefs(FS_ori_1.mod.red2, FS_ori_1.mod.red) # compares estimated coefficients and their standard errors
#It doesn't look like they change much

Anova(FS_ori_1.mod.red)

effect_plot(FS_ori_1.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'ENN ', y.label = 'Standardized Effect Sizes (SES)')

FS_ori_1.mod.null <- glm(SES_ori_1 ~ 1, family = gaussian, data = FS_comps)

# model comparison techniques
anova(FS_ori_1.mod.full, FS_ori_1.mod.red, FS_ori_1.mod.null)
AICctab(FS_ori_1.mod.full, FS_ori_1.mod.red, FS_ori_1.mod.null)


KT_ori_1.mod.full <- glm(SES_ori_1 ~ trmt + pland + enn, family = gaussian, data = KT_comps)
summary(KT_ori_1.mod.full)
step(KT_ori_1.mod.full)


KT_ori_1.mod.null <- glm(SES_ori_0 ~ 1, family = gaussian, data = KT_comps)
summary(KT_ori_1.mod.null)
qqnorm(resid(KT_ori_1.mod.null))
qqline(resid(KT_ori_1.mod.null))
plot(simulateResiduals(KT_ori_1.mod.null))
densityPlot(rstudent(KT_ori_1.mod.null)) # check density estimate of the distribution of residuals
outlierTest(KT_ori_1.mod.null)
influenceIndexPlot(KT_ori_1.mod.null, vars = c("Cook"), id = list(n = 3))


# model comparison techniques
anova(KT_ori_1.mod.full, KT_ori_1.mod.null)
AICctab(KT_ori_1.mod.full,  KT_ori_1.mod.null)

VL_ori_1.mod.full <- glm(SES_ori_1 ~ trmt + pland + enn, family = gaussian, data = VL_comps)
summary(VL_ori_1.mod.full)
step(VL_ori_1.mod.full)


VL_ori_1.mod.red <- glm(SES_ori_1 ~ enn, family = gaussian, data = VL_comps)
summary(VL_ori_1.mod.red)
qqnorm(resid(VL_ori_1.mod.red))
qqline(resid(VL_ori_1.mod.red))
plot(simulateResiduals(VL_ori_1.mod.red))
#quantile deviations detected
densityPlot(rstudent(VL_ori_1.mod.red)) # check density estimate of the distribution of residuals
outlierTest(VL_ori_1.mod.red)
influenceIndexPlot(VL_ori_1.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(VL_ori_1.mod.red)

effect_plot(VL_ori_1.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Largest Patch Index (LPI)', y.label = 'Standardized Effect Sizes (SES)')

VL_ori_1.mod.null <- glm(SES_ori_1 ~ 1, family = gaussian, data = VL_comps)

# model comparison techniques
anova(VL_ori_1.mod.full, VL_ori_1.mod.red, VL_ori_1.mod.null)
AICctab(VL_ori_1.mod.full, VL_ori_1.mod.red, VL_ori_1.mod.null)

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
with(FS_comps, bartlett.test(SES_bsor ~ trmt))
with(KT_comps, bartlett.test(SES_bsor ~ trmt))
with(VL_comps, bartlett.test(SES_bsor ~ trmt))
with(FS_comps, ad.test(SES_bsor))
with(KT_comps, ad.test(SES_bsor))
with(VL_comps, ad.test(SES_bsor))


FS_bsor.mod.full <- glm(SES_bsor ~ trmt + pland + enn, family = gaussian, data = FS_comps)
summary(FS_bsor.mod.full)
step(FS_bsor.mod.full)

FS_bsor.mod.null <- glm(SES_bsor ~ 1, family = gaussian, data = FS_comps)
summary(FS_bsor.mod.null)
qqnorm(resid(FS_bsor.mod.null))
qqline(resid(FS_bsor.mod.null))
plot(simulateResiduals(FS_bsor.mod.null))
densityPlot(rstudent(FS_bsor.mod.null)) # check density estimate of the distribution of residuals
outlierTest(FS_bsor.mod.null)
influenceIndexPlot(FS_soc_4.mod.null, vars = c("Cook"), id = list(n = 3))


# model comparison techniques
anova(FS_bsor.mod.full, FS_bsor.mod.null)
AICctab(FS_bsor.mod.full, FS_bsor.mod.null)


KT_bsor.mod.full <- glm(SES_bsor ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_bsor.mod.full)
step(KT_bsor.mod.full)

KT_bsor.mod.red <- glm(SES_bsor ~ enn, family = gaussian, data = KT_comps)
summary(KT_bsor.mod.red)
qqnorm(resid(KT_bsor.mod.red))
qqline(resid(KT_bsor.mod.red))
plot(simulateResiduals(KT_bsor.mod.red))
#Quantile deviations detected
densityPlot(rstudent(KT_bsor.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_bsor.mod.red)
influenceIndexPlot(KT_bsor.mod.red, vars = c("Cook"), id = list(n = 3))

KT_bsor.mod.null <- glm(SES_bsor ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_bsor.mod.full, KT_bsor.mod.red, KT_bsor.mod.null)
AICctab(KT_bsor.mod.full, KT_bsor.mod.red, KT_bsor.mod.null)

#Finding effect values of model
Anova(KT_bsor.mod.red)

effect_plot(KT_bsor.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Greenspace Isolation (ENN)', y.label = 'Standardized Effect Sizes (SES)')


VL_bsor.mod.full <- glm(SES_bsor ~ trmt + pland + enn, family = gaussian, data = VL_comps)
summary(VL_bsor.mod.full)
step(VL_bsor.mod.full)

VL_bsor.mod.red <- glm(SES_bsor ~ enn, family = gaussian, data = VL_comps)
summary(VL_bsor.mod.red)
qqnorm(resid(VL_bsor.mod.red))
qqline(resid(VL_bsor.mod.red))
plot(simulateResiduals(VL_bsor.mod.red))
densityPlot(rstudent(VL_bsor.mod.red)) # check density estimate of the distribution of residuals
outlierTest(VL_bsor.mod.red)
influenceIndexPlot(VL_bsor.mod.red, vars = c("Cook"), id = list(n = 3))

VL_bsor.mod.null <- glm(SES_bsor ~ 1, family = gaussian, data = VL_comps)

# model comparison techniques
anova(VL_bsor.mod.full, VL_bsor.mod.red, VL_bsor.mod.null)
AICctab(VL_bsor.mod.full, VL_bsor.mod.red, VL_bsor.mod.null)

#Finding effect values of model
Anova(VL_bsor.mod.red)

effect_plot(VL_bsor.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Greenspace Isolation (ENN)', y.label = 'Standardized Effect Sizes (SES)')


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
bsim.t1 <- wilcox.test(T1$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsim.t1
bsim.t8 <- wilcox.test(T8$SES_bsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsim.t8


## compare among treatments and landscape variables
dotchart(SES$SES_bsim, group = SES$trmt, pch = 19)
with(FS_comps, bartlett.test(SES_bsim ~ trmt))
with(KT_comps, bartlett.test(SES_bsim ~ trmt))
with(VL_comps, bartlett.test(SES_bsim ~ trmt))
with(FS_comps, ad.test(SES_bsim))
with(KT_comps, ad.test(SES_bsim))
with(VL_comps, ad.test(SES_bsim))


FS_bsim.mod.full <- glm(SES_bsim ~ trmt + pland  + enn, family = gaussian, data = FS_comps)
summary(FS_bsim.mod.full)
step(FS_bsim.mod.full)

FS_bsim.mod.red <- glm(SES_bsim ~ trmt+enn, family = gaussian, data = FS_comps)
summary(FS_bsim.mod.red)
qqnorm(resid(FS_bsim.mod.red))
qqline(resid(FS_bsim.mod.red))
plot(simulateResiduals(FS_bsim.mod.red))
densityPlot(rstudent(FS_bsim.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_bsim.mod.red)
influenceIndexPlot(FS_bsim.mod.red, vars = c("Cook"), id = list(n = 3))

FS_bsim.mod.null <- glm(SES_bsim ~ 1, family = gaussian, data = FS_comps)

# model comparison techniques
anova(FS_bsim.mod.full, FS_bsim.mod.red, FS_bsim.mod.null)
AICctab(FS_bsim.mod.full, FS_bsim.mod.red, FS_bsim.mod.null)
#null was best model

KT_bsim.mod.full <- glm(SES_bsim ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_bsim.mod.full)
step(KT_bsim.mod.full)

KT_bsim.mod.red <- glm(SES_bsim ~ pland, family = gaussian, data = KT_comps)
summary(KT_bsim.mod.red)
qqnorm(resid(KT_bsim.mod.red))
qqline(resid(KT_bsim.mod.red))
plot(simulateResiduals(KT_bsim.mod.red))
#Quantile deviations detected
densityPlot(rstudent(KT_bsim.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_bsim.mod.red)
influenceIndexPlot(KT_bsim.mod.red, vars = c("Cook"), id = list(n = 3))

KT_bsim.mod.null <- glm(SES_bsim ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_bsim.mod.full, KT_bsim.mod.red, KT_bsim.mod.null)
AICctab(KT_bsim.mod.full, KT_bsim.mod.red, KT_bsim.mod.null)

#null is best model

VL_bsim.mod.full <- glm(SES_bsim ~ trmt + pland + enn, family = gaussian, data = VL_comps)
summary(VL_bsim.mod.full)
step(VL_bsim.mod.full)

VL_bsim.mod.red <- glm(SES_bsim ~ trmt + enn, family = gaussian, data = VL_comps)
summary(VL_bsim.mod.red)
qqnorm(resid(VL_bsim.mod.red))
qqline(resid(VL_bsim.mod.red))
plot(simulateResiduals(VL_bsim.mod.red))
#no sig problems detected
densityPlot(rstudent(VL_bsim.mod.red)) # check density estimate of the distribution of residuals
outlierTest(VL_bsim.mod.red)
influenceIndexPlot(VL_bsim.mod.red, vars = c("Cook"), id = list(n = 3))

VL_bsim.mod.null <- glm(SES_bsim ~ 1, family = gaussian, data = VL_comps)

# model comparison techniques
anova(VL_bsim.mod.full, VL_bsim.mod.red, VL_bsim.mod.null)
AICctab(VL_bsim.mod.full, VL_bsim.mod.red, VL_bsim.mod.null)

#null was best model

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
bsne.t1 <- wilcox.test(T1$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsne.t1
bsne.t8 <- wilcox.test(T8$SES_bsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
bsne.t8


## compare among treatments and landscape variables
dotchart(SES$SES_bsne, group = SES$trmt, pch = 19)
with(FS_comps, bartlett.test(SES_bsne ~ trmt))
with(KT_comps, bartlett.test(SES_bsne ~ trmt))
with(VL_comps, bartlett.test(SES_bsne ~ trmt))
with(FS_comps, ad.test(SES_bsne))
with(KT_comps, ad.test(SES_bsne))
with(VL_comps, ad.test(SES_bsne))
#VL fails Anderson Darling normality test. KW test instead

FS_bsne.mod.full <- glm(SES_bsne ~ trmt + pland  + enn, family = gaussian, data = FS_comps)
summary(FS_bsne.mod.full)
step(FS_bsne.mod.full)

FS_bsne.mod.red <- glm(SES_bsne ~ trmt, family = gaussian, data = FS_comps)
summary(FS_bsne.mod.red)
qqnorm(resid(FS_bsne.mod.red))
qqline(resid(FS_bsne.mod.red))
plot(simulateResiduals(FS_bsne.mod.red))
densityPlot(rstudent(FS_bsne.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_bsne.mod.red)
influenceIndexPlot(FS_bsne.mod.red, vars = c("Cook"), id = list(n = 3))

FS_bsne.mod.null <- glm(SES_bsne ~ 1, family = gaussian, data = FS_comps)

# model comparison techniques
anova(FS_bsne.mod.full, FS_bsne.mod.red, FS_bsne.mod.null)
AICctab(FS_bsne.mod.full, FS_bsne.mod.red, FS_bsne.mod.null)

#Null was best model fit

KT_bsne.mod.full <- glm(SES_bsne ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_bsne.mod.full)
step(KT_bsne.mod.full)

KT_bsne.mod.red <- glm(SES_bsne ~ pland, family = gaussian, data = KT_comps)
summary(KT_bsne.mod.red)
qqnorm(resid(KT_bsne.mod.red))
qqline(resid(KT_bsne.mod.red))
plot(simulateResiduals(KT_bsne.mod.red))
#quantile deviations again
densityPlot(rstudent(KT_bsne.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_bsne.mod.red)
influenceIndexPlot(KT_bsne.mod.red, vars = c("Cook"), id = list(n = 3))

KT_bsne.mod.null <- glm(SES_bsne ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_bsne.mod.full, KT_bsne.mod.red, KT_bsne.mod.null)
AICctab(KT_bsne.mod.full, KT_bsne.mod.red, KT_bsne.mod.null)

#Null was best model fit

#kruskal wallace for VL
kruskal.test(SES_bsne ~ trmt, data = VL_comps)


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
with(FS_comps, bartlett.test(SES_falpha ~ trmt))
with(KT_comps, bartlett.test(SES_falpha ~ trmt))
with(VL_comps, bartlett.test(SES_falpha ~ trmt))
with(FS_comps, ad.test(SES_falpha))
with(KT_comps, ad.test(SES_falpha))
with(VL_comps, ad.test(SES_falpha))


FS_falpha.mod.full <- glm(SES_falpha ~ trmt + pland  + enn, family = gaussian, data = FS_comps)
summary(FS_falpha.mod.full)
step(FS_falpha.mod.full)

FS_falpha.mod.null <- glm(SES_falpha ~ 1, family = gaussian, data = FS_comps)
summary(FS_falpha.mod.null)
qqnorm(resid(FS_falpha.mod.null))
qqline(resid(FS_falpha.mod.null))
plot(simulateResiduals(FS_falpha.mod.null))
densityPlot(rstudent(FS_falpha.mod.null)) # check density estimate of the distribution of residuals
outlierTest(FS_falpha.mod.null)
influenceIndexPlot(FS_falpha.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(FS_falpha.mod.full, FS_falpha.mod.null)
AICtab(FS_falpha.mod.full, FS_falpha.mod.null)


KT_falpha.mod.full <- glm(SES_falpha ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_falpha.mod.full)
step(KT_falpha.mod.full)

KT_falpha.mod.red <- glm(SES_falpha ~ trmt, family = gaussian, data = KT_comps)
summary(KT_falpha.mod.red)
qqnorm(resid(KT_falpha.mod.red))
qqline(resid(KT_falpha.mod.red))
plot(simulateResiduals(KT_falpha.mod.red))
densityPlot(rstudent(KT_falpha.mod.red)) # check density estimate of the distribution of residuals
outlierTest(KT_falpha.mod.red)
influenceIndexPlot(KT_falpha.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(KT_falpha.mod.red)
emmeans(KT_falpha.mod.red, pairwise ~ trmt)

effect_plot(KT_falpha.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = '', y.label = 'Standardized Effect Sizes (SES)')

KT_falpha.mod.null <- glm(SES_falpha ~ 1, family = gaussian, data = KT_comps)

# model comparison techniques
anova(KT_falpha.mod.full, KT_falpha.mod.red, KT_falpha.mod.null)
AICtab(KT_falpha.mod.full, KT_falpha.mod.red, KT_falpha.mod.null)


VL_falpha.mod.full <- glm(SES_falpha ~ trmt + pland  + enn, family = gaussian, data = VL_comps)
summary(VL_falpha.mod.full)
step(VL_falpha.mod.full)

VL_falpha.mod.red <- glm(SES_falpha ~ pland + enn, family = gaussian, data = VL_comps)
summary(VL_falpha.mod.red)
qqnorm(resid(VL_falpha.mod.red))
qqline(resid(VL_falpha.mod.red))
plot(simulateResiduals(VL_falpha.mod.red))
densityPlot(rstudent(VL_falpha.mod.red)) # check density estimate of the distribution of residuals
outlierTest(VL_falpha.mod.red)
#There's an outlier. 
influenceIndexPlot(VL_falpha.mod.red, vars = c("Cook"), id = list(n = 3))
#None of the points are significantly influencing the model though.
#I'm going to skip removing the outlier, because it's not really influencing the model

Anova(VL_falpha.mod.red)

effect_plot(VL_falpha.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = '(PLAND)', y.label = 'SES')
effect_plot(VL_falpha.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'ENN', y.label = 'SES')


VL_falpha.mod.null <- glm(SES_falpha ~ 1, family = gaussian, data = VL_comps)

# model comparison techniques
anova(VL_falpha.mod.full, VL_falpha.mod.red, VL_falpha.mod.null)
AICtab(VL_falpha.mod.full, VL_falpha.mod.red, VL_falpha.mod.null)


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
with(FS_comps, bartlett.test(SES_fbsor ~ trmt))
with(KT_comps, bartlett.test(SES_fbsor ~ trmt))
with(VL_comps, bartlett.test(SES_fbsor ~ trmt))
with(FS_comps, ad.test(SES_fbsor))
with(KT_comps, ad.test(SES_fbsor))
with(VL_comps, ad.test(SES_fbsor))


FS_fbsor.mod.full <- glm(SES_fbsor ~ trmt + pland + enn, family = gaussian, data = FS_comps)
summary(FS_fbsor.mod.full)
step(FS_fbsor.mod.full)

FS_fbsor.mod.red<-glm(SES_fbsor ~ trmt, family = gaussian, data = FS_comps)
FS_fbsor.mod.null <- glm(SES_fbsor ~ 1, family = gaussian, data = FS_comps)
summary(FS_fbsor.mod.red)
qqnorm(resid(FS_fbsor.mod.red))
qqline(resid(FS_fbsor.mod.red))
plot(simulateResiduals(FS_fbsor.mod.red))
densityPlot(rstudent(FS_fbsor.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_fbsor.mod.red)

influenceIndexPlot(FS_fbsor.mod.red, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(FS_fbsor.mod.full,FS_fbsor.mod.red, FS_fbsor.mod.null)
AICtab(FS_fbsor.mod.full, FS_fbsor.mod.red, FS_fbsor.mod.null)

Anova(FS_fbsor.mod.red)

png("Figure trmt fbsor.png", width = 1500, height = 1000, pointsize = 20)
effect_plot(FS_fbsor.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'Treatment Comparison', y.label = 'Standardized Effect Sizes (SES)')
dev.off()


KT_fbsor.mod.full <- glm(SES_fbsor ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_fbsor.mod.full)
step(KT_fbsor.mod.full)

KT_fbsor.mod.null <- glm(SES_fbsor ~ 1, family = gaussian, data = KT_comps)
summary(KT_fbsor.mod.null)
qqnorm(resid(KT_fbsor.mod.null))
qqline(resid(KT_fbsor.mod.null))
plot(simulateResiduals(KT_fbsor.mod.null))
densityPlot(rstudent(KT_fbsor.mod.null)) # check density estimate of the distribution of residuals
outlierTest(KT_fbsor.mod.null)
influenceIndexPlot(KT_fbsor.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(KT_fbsor.mod.full, KT_fbsor.mod.null)
AICtab(KT_fbsor.mod.full, KT_fbsor.mod.null)


VL_fbsor.mod.full <- glm(SES_fbsor ~ trmt + pland  + enn, family = gaussian, data = VL_comps)
summary(VL_fbsor.mod.full)
step(VL_fbsor.mod.full)

VL_fbsor.mod.null <- glm(SES_fbsor ~ 1, family = gaussian, data = VL_comps)
summary(VL_fbsor.mod.null)
qqnorm(resid(VL_fbsor.mod.null))
qqline(resid(VL_fbsor.mod.null))
plot(simulateResiduals(VL_fbsor.mod.null))
densityPlot(rstudent(VL_fbsor.mod.null)) # check density estimate of the distribution of residuals
outlierTest(VL_fbsor.mod.null)
influenceIndexPlot(VL_fbsor.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(VL_fbsor.mod.full, VL_fbsor.mod.null)
AICtab(VL_fbsor.mod.full, VL_fbsor.mod.null)


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
fbsim.t1 <- wilcox.test(T1$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsim.t1
fbsim.t8 <- wilcox.test(T8$SES_fbsim, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsim.t8


## compare among treatments and landscape variables
dotchart(SES$SES_fbsim, group = SES$trmt, pch = 19)
with(FS_comps, bartlett.test(SES_fbsim ~ trmt))
with(KT_comps, bartlett.test(SES_fbsim ~ trmt))
with(VL_comps, bartlett.test(SES_fbsim ~ trmt))
with(FS_comps, ad.test(SES_fbsim))
with(KT_comps, ad.test(SES_fbsim))
with(VL_comps, ad.test(SES_fbsim))
 
FS_fbsim.mod.full <- glm(SES_fbsim ~ trmt + pland + enn, family = gaussian, data = FS_comps)
summary(FS_fbsim.mod.full)
step(FS_fbsim.mod.full)

FS_fbsim.mod.red <- glm(SES_fbsim ~ pland , family = gaussian, data = FS_comps)
summary(FS_fbsim.mod.red)
qqnorm(resid(FS_fbsim.mod.red))
qqline(resid(FS_fbsim.mod.red))
plot(simulateResiduals(FS_fbsim.mod.red))
densityPlot(rstudent(FS_fbsim.mod.red)) # check density estimate of the distribution of residuals
outlierTest(FS_fbsim.mod.red)
influenceIndexPlot(FS_fbsim.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(FS_fbsim.mod.red)
effect_plot(FS_fbsim.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = '(PLAND)', y.label = 'SES')


FS_fbsim.mod.null <- glm(SES_fbsim ~ 1, family = gaussian, data = FS_comps)

# model comparison techniques
anova(FS_fbsim.mod.full, FS_fbsim.mod.red, FS_fbsim.mod.null)
AICtab(FS_fbsim.mod.full, FS_fbsim.mod.red, FS_fbsim.mod.null)


KT_fbsim.mod.full <- glm(SES_fbsim ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_fbsim.mod.full)
step(KT_fbsim.mod.full)

KT_fbsim.mod.null <- glm(SES_fbsim ~ 1, family = gaussian, data = KT_comps)
summary(KT_fbsim.mod.null)
qqnorm(resid(KT_fbsim.mod.null))
qqline(resid(KT_fbsim.mod.null))
plot(simulateResiduals(KT_fbsim.mod.null))
densityPlot(rstudent(KT_fbsim.mod.null)) # check density estimate of the distribution of residuals
outlierTest(KT_fbsim.mod.null)
influenceIndexPlot(KT_fbsim.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(KT_fbsim.mod.full, KT_fbsim.mod.null)
AICtab(KT_fbsim.mod.full, KT_fbsim.mod.null)


VL_fbsim.mod.full <- glm(SES_fbsim ~ trmt + pland  + enn, family = gaussian, data = VL_comps)
summary(VL_fbsim.mod.full)
step(VL_fbsim.mod.full)

VL_fbsim.mod.red <- glm(SES_fbsim ~ pland + enn, family = gaussian, data = VL_comps)
summary(VL_fbsim.mod.red)
qqnorm(resid(VL_fbsim.mod.red))
qqline(resid(VL_fbsim.mod.red))
plot(simulateResiduals(VL_fbsim.mod.red))
densityPlot(rstudent(VL_fbsim.mod.red)) # check density estimate of the distribution of residuals
outlierTest(VL_fbsim.mod.red)
influenceIndexPlot(VL_fbsim.mod.red, vars = c("Cook"), id = list(n = 3))

Anova(VL_fbsim.mod.red)
effect_plot(VL_fbsim.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = '(PLAND)', y.label = '')
effect_plot(VL_fbsim.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = '(ENN)', y.label = ' ')


VL_fbsim.mod.null <- glm(SES_fbsim ~ 1, family = gaussian, data = VL_comps)

# model comparison techniques
anova(VL_fbsim.mod.full, VL_fbsim.mod.red, VL_fbsim.mod.null)
AICtab(VL_fbsim.mod.full, VL_fbsim.mod.red, VL_fbsim.mod.null)



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
fbsne.t1 <- wilcox.test(T1$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsne.t1
fbsne.t8 <- wilcox.test(T8$SES_fbsne, y = NULL, mu = 0, alternative = c("two.sided"), conf.int = TRUE)
fbsne.t8


## compare among treatments and landscape variables
dotchart(SES$SES_fbsne, group = SES$trmt, pch = 19)
with(FS_comps, bartlett.test(SES_fbsne ~ trmt))
with(KT_comps, bartlett.test(SES_fbsne ~ trmt))
with(VL_comps, bartlett.test(SES_fbsne ~ trmt))
with(FS_comps, ad.test(SES_fbsne))
#Failed normality test. I'm not seeing any major outliers?
with(KT_comps, ad.test(SES_fbsne))
with(VL_comps, ad.test(SES_fbsne))

#Skipping FS because of failed normality test

KT_fbsne.mod.full <- glm(SES_fbsne ~ trmt + pland  + enn, family = gaussian, data = KT_comps)
summary(KT_fbsne.mod.full)
step(KT_fbsne.mod.full)

KT_fbsne.mod.null <- glm(SES_fbsne ~ 1, family = gaussian, data = KT_comps)
summary(KT_fbsne.mod.null)
qqnorm(resid(KT_fbsne.mod.null))
qqline(resid(KT_fbsne.mod.null))
plot(simulateResiduals(KT_fbsne.mod.null))
densityPlot(rstudent(KT_fbsne.mod.null)) # check density estimate of the distribution of residuals
outlierTest(KT_fbsne.mod.null)
influenceIndexPlot(KT_fbsne.mod.null, vars = c("Cook"), id = list(n = 3))

# model comparison techniques
anova(KT_fbsne.mod.full, KT_fbsne.mod.null)
AICctab(KT_fbsne.mod.full, KT_fbsne.mod.null)


VL_fbsne.mod.full <- glm(SES_fbsne ~ trmt + pland + enn, family = gaussian, data = VL_comps)
summary(VL_fbsne.mod.full)
step(VL_fbsne.mod.full)

VL_fbsne.mod.null <- glm(SES_fbsne ~ 1, family = gaussian, data = VL_comps)
summary(VL_fbsne.mod.null)
qqnorm(resid(VL_fbsne.mod.null))
qqline(resid(VL_fbsne.mod.null))
plot(simulateResiduals(VL_fbsne.mod.null))
densityPlot(rstudent(VL_fbsne.mod.null)) # check density estimate of the distribution of residuals
outlierTest(VL_fbsne.mod.null)
influenceIndexPlot(VL_fbsne.mod.null, vars = c("Cook"), id = list(n = 3))


# model comparison techniques
anova(VL_fbsne.mod.full, VL_fbsne.mod.null)
AICtab(VL_fbsne.mod.full, VL_fbsne.mod.null)



#Graph making----
library(ggplot2)
#install.packages("ggthemes")
library(ggthemes)


boxplot(SES_bsor ~ trmt, data = FS_comps, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, cex.lab = 2, cex.axis = 1.7)
stripchart(ses ~ fbeta, data = SES_TBeta, col = viridis(3),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(52, 3.4, "A", pos = 4, font = 2, cex = 2.6)


boxplot(ses ~ fbeta, data = FS_comps, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, cex.lab = 2, cex.axis = 1.7)
stripchart(ses ~ fbeta, data = SES_FBeta, col = viridis(3),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(5.5, 3.4, "B", pos = 4, font = 2, cex = 2.6)

dev.off()

#text(1700, 17.5, "A", pos = 4, font = 2, cex = 3)

