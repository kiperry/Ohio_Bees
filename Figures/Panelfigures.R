#Making Panel figures 
#Edited by Caralee Shepard 8.8.24
#edited again by C.S. 1.15.25 (just to change figure colors)
#editing again 1.17.25 following changes with updated species list

#Step 1- importing data----
##1.1 Pham data----
a.mp <- read.csv("urbanpool.MP.26.csv", row.names=1)

aO.mp<-a.mp

a.mp<-a.mp[1:136]
rowSums(a.mp)



#import the SES data
SES.mp <- read.csv("Urban to Local_Nulls_MP/SES_Local_MP.csv", row.names = 1)

#Including treatment
a1.mp<-aO.mp


SES.mp$trmt <- a1.mp$trmt
str(SES.mp)
SES.mp

#Adding the random neighborhood variable
SES.mp$neighd <- c("C", "DS", "F", "G", "SV","TM","SV","TM", "H", "SV", "BE", "C","DS","F","G")


## pull out data for each treatment

VL2019<- SES.mp[which(SES.mp$trmt== "VL"),]
str(VL2019)
UP2019<- SES.mp[which(SES.mp$trmt== "UP"),]
str(UP2019)

##1.2 Turo data ----
a.kt <- read.csv("urbanpool.KT.26.csv", row.names=1)

aO.kt<-a.kt
a.kt<-a.kt[1:136]

#import the SES data
SES.kt <- read.csv("Urban to Local_Nulls_KT/SES_Local_KT.csv", row.names = 1)
SES.kt

#Including treatment

SES.kt$trmt <- aO.kt$trmt
str(SES.kt)

#Adding the random neighborhood variable
SES.kt$neighd <- c("BU", "CE", "DS", "FA", "GL","HO","SV","TR","BU", "CE", "DS", "FA", "GL","HO","SV","TR","BU", "CE", "DS", "FA", "GL","HO","SV","TR")


## pull out data for each treatment
T1 <- SES.kt[which(SES.kt$trmt == "T1"),]
str(T1)
Prairie <- SES.kt[which(SES.kt$trmt == "Prairie"),]
str(Prairie)

##1.3 Sivakoff data----
a.fs <- read.csv("urbanpool.FS.26.csv", row.names=1)

aO.fs<-a.fs
a.fs<-a.fs[1:136]

#import the SES data
SES.fs <- read.csv("Urban to Local_Nulls_FS/SES_Local_FS.csv", row.names = 1)
SES.fs
#Including treatment

SES.fs$trmt <- aO.fs$trmt
str(SES.fs)

## pull out data for each treatment
farm <- SES.fs[which(SES.fs$trmt == "Farm"),]
str(farm)

control <- SES.fs[which(SES.fs$trmt == "Control"),]
str(control)





#Figures----
##Loading programs needed----
if (!suppressWarnings(require(viridis))) install.packages("viridis")
citation("viridis")

if (!suppressWarnings(require(reshape2))) install.packages("reshape2")
citation("reshape2")

##Taxonomic and functional diversity panel----
###Option 1----
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


SES_ALLdivt1 <- as.data.frame(rbind( T1$SES_fbsim, T1$SES_fbsne, T1$SES_fbsor, T1$SES_bsim, T1$SES_bsne, T1$SES_bsor, T1$SES_falpha))
str(SES_ALLdivt1)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdivt1 <- data.frame(t(SES_ALLdivt1))
str(SES_ALLdivt1)
colnames(SES_ALLdivt1) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity","Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity", "Funct. Alpha" )


SES_ALLdivt1 <- melt(SES_ALLdivt1)
colnames(SES_ALLdivt1) <- c("diversity","ses")



SES_ALLdivprai <- as.data.frame(rbind( Prairie$SES_fbsim, Prairie$SES_fbsne, Prairie$SES_fbsor, Prairie$SES_bsim, Prairie$SES_bsne, Prairie$SES_bsor, Prairie$SES_falpha))
str(SES_ALLdivprai)
SES_ALLdivprai <- data.frame(t(SES_ALLdivprai))
str(SES_ALLdivprai)
colnames(SES_ALLdivprai) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity","Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity", "Funct. Alpha")

SES_ALLdivprai <- melt(SES_ALLdivprai)
colnames(SES_ALLdivprai) <- c("diversity","ses")





#Now I should be able to plot them in a nice graph
png("Figures/paneltaxandfunct.option1.png", width = 2500, height = 2500, pointsize = 20)
par(mfrow=c(3,2)) # indicates 3 rows, 2 columns
par(mar = c(5,11,4,2)) # sets the margins around the figure

boxplot(ses ~ diversity, data = SES_ALLdivcontrol, col = magma(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-27,25),cex.lab = 2, cex.axis = 1.45)
stripchart(ses ~ diversity, data = SES_ALLdivcontrol, col = viridis(6),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(20, 7, "A", pos = 4, font = 2, cex = 2.6)


boxplot(ses ~ diversity, data = SES_ALLdivfarm, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-27,22.5))
stripchart(ses ~ diversity, data = SES_ALLdivfarm, col = viridis(6),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(20, 7, "B", pos = 4, font = 2, cex = 2.6)


boxplot(ses ~ diversity, data = SES_ALLdivt1, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-27,25))
stripchart(ses ~ diversity, data = SES_ALLdivt1, col = viridis(6),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(20, 7, "C", pos = 4, font = 2, cex = 2.6)


boxplot(ses ~ diversity, data = SES_ALLdivprai, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-27,22.5))
stripchart(ses ~ diversity, data = SES_ALLdivprai, col = viridis(6),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(20, 7, "D", pos = 4, font = 2, cex = 2.6)



boxplot(ses ~ diversity, data = SES_ALLdivvl, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-27,22.5))
stripchart(ses ~ diversity, data = SES_ALLdivvl, col = viridis(6),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(20, 7, "E", pos = 4, font = 2, cex = 2.6)


boxplot(ses ~ diversity, data = SES_ALLdivup, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-27,22.5))
stripchart(ses ~ diversity, data = SES_ALLdivup, col = viridis(6),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(20, 7, "F", pos = 4, font = 2, cex = 2.6)

dev.off()
###Option 2----
#step one- create a dataframe for all of the tax and funct diversity metrics I want to graph
SES_ALLdivvl.tax <- as.data.frame(rbind(VL2019$SES_bsim, VL2019$SES_bsne, VL2019$SES_bsor))
str(SES_ALLdivvl.tax)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdivvl.tax <- data.frame(t(SES_ALLdivvl.tax))
str(SES_ALLdivvl.tax)
colnames(SES_ALLdivvl.tax) <- c("Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity" )


SES_ALLdivvl.tax <- melt(SES_ALLdivvl.tax)
colnames(SES_ALLdivvl.tax) <- c("diversity","ses")


SES_ALLdivvl.func <- as.data.frame(rbind( VL2019$SES_fbsim, VL2019$SES_fbsne, VL2019$SES_fbsor, VL2019$SES_falpha))
str(SES_ALLdivvl.func)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdivvl.func  <- data.frame(t(SES_ALLdivvl.func ))
str(SES_ALLdivvl.func )
colnames(SES_ALLdivvl.func ) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity", "Funct. Alpha" )


SES_ALLdivvl.func  <- melt(SES_ALLdivvl.func )
colnames(SES_ALLdivvl.func) <- c("diversity","ses")



SES_ALLdivup.tax <- as.data.frame(rbind( UP2019$SES_bsim, UP2019$SES_bsne, UP2019$SES_bsor))
str(SES_ALLdivup.tax)
SES_ALLdivup.tax <- data.frame(t(SES_ALLdivup.tax))
str(SES_ALLdivup.tax)
colnames(SES_ALLdivup.tax) <- c("Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity" )

SES_ALLdivup.tax <- melt(SES_ALLdivup.tax)
colnames(SES_ALLdivup.tax) <- c("diversity","ses")

SES_ALLdivup.func <- as.data.frame(rbind( UP2019$SES_fbsim, UP2019$SES_fbsne, UP2019$SES_fbsor, UP2019$SES_falpha))
str(SES_ALLdivup.func)
SES_ALLdivup.func <- data.frame(t(SES_ALLdivup.func))
str(SES_ALLdivup.func)
colnames(SES_ALLdivup.func) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity", "Funct. Alpha" )

SES_ALLdivup.func <- melt(SES_ALLdivup.func)
colnames(SES_ALLdivup.func) <- c("diversity","ses")

####sivakoff dataframes-----

SES_ALLdivcontrol.tax <- as.data.frame(rbind( control$SES_bsim, control$SES_bsne, control$SES_bsor))
str(SES_ALLdivcontrol.tax)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdivcontrol.tax <- data.frame(t(SES_ALLdivcontrol.tax))
str(SES_ALLdivcontrol.tax)
colnames(SES_ALLdivcontrol.tax) <- c("Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity" )


SES_ALLdivcontrol.tax <- melt(SES_ALLdivcontrol.tax)
colnames(SES_ALLdivcontrol.tax) <- c("diversity","ses")

SES_ALLdivcontrol.func <- as.data.frame(rbind( control$SES_fbsim, control$SES_fbsne, control$SES_fbsor, control$SES_falpha))
str(SES_ALLdivcontrol.func)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdivcontrol.func <- data.frame(t(SES_ALLdivcontrol.func))
str(SES_ALLdivcontrol.func)
colnames(SES_ALLdivcontrol.func) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity", "Funct. Alpha" )


SES_ALLdivcontrol.func <- melt(SES_ALLdivcontrol.func)
colnames(SES_ALLdivcontrol.func) <- c("diversity","ses")


SES_ALLdivfarm.tax <- as.data.frame(rbind( farm$SES_bsim, farm$SES_bsne, farm$SES_bsor))
str(SES_ALLdivfarm.tax)
SES_ALLdivfarm.tax <- data.frame(t(SES_ALLdivfarm.tax))
str(SES_ALLdivfarm.tax)
colnames(SES_ALLdivfarm.tax) <- c("Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity" )

SES_ALLdivfarm.tax <- melt(SES_ALLdivfarm.tax)
colnames(SES_ALLdivfarm.tax) <- c("diversity","ses")

SES_ALLdivfarm.func <- as.data.frame(rbind( farm$SES_fbsim, farm$SES_fbsne, farm$SES_fbsor, farm$SES_falpha))
str(SES_ALLdivfarm.func)
SES_ALLdivfarm.func <- data.frame(t(SES_ALLdivfarm.func))
str(SES_ALLdivfarm.func)
colnames(SES_ALLdivfarm.func) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity", "Funct. Alpha" )

SES_ALLdivfarm.func <- melt(SES_ALLdivfarm.func)
colnames(SES_ALLdivfarm.func) <- c("diversity","ses")

####Turodataframes----
SES_ALLdivt1.tax <- as.data.frame(rbind( T1$SES_bsim, T1$SES_bsne, T1$SES_bsor))
str(SES_ALLdivt1.tax)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdivt1.tax <- data.frame(t(SES_ALLdivt1.tax))
str(SES_ALLdivt1.tax)
colnames(SES_ALLdivt1.tax) <- c("Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity" )


SES_ALLdivt1.tax <- melt(SES_ALLdivt1.tax)
colnames(SES_ALLdivt1.tax) <- c("diversity","ses")

SES_ALLdivt1.func <- as.data.frame(rbind( T1$SES_fbsim, T1$SES_fbsne, T1$SES_fbsor, T1$SES_falpha))
str(SES_ALLdivt1.func)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdivt1.func <- data.frame(t(SES_ALLdivt1.func))
str(SES_ALLdivt1.func)
colnames(SES_ALLdivt1.func) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity", "Funct. Alpha" )


SES_ALLdivt1.func <- melt(SES_ALLdivt1.func)
colnames(SES_ALLdivt1.func) <- c("diversity","ses")


SES_ALLdivprai.tax <- as.data.frame(rbind( Prairie$SES_bsim, Prairie$SES_bsne, Prairie$SES_bsor))
str(SES_ALLdivprai.tax)
SES_ALLdivprai.tax <- data.frame(t(SES_ALLdivprai.tax))
str(SES_ALLdivprai.tax)
colnames(SES_ALLdivprai.tax) <- c("Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity" )

SES_ALLdivprai.tax <- melt(SES_ALLdivprai.tax)
colnames(SES_ALLdivprai.tax) <- c("diversity","ses")

SES_ALLdivprai.func <- as.data.frame(rbind( Prairie$SES_fbsim, Prairie$SES_fbsne, Prairie$SES_fbsor, Prairie$SES_falpha))
str(SES_ALLdivprai.func)
SES_ALLdivprai.func <- data.frame(t(SES_ALLdivprai.func))
str(SES_ALLdivprai.func)
colnames(SES_ALLdivprai.func) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity", "Funct. Alpha" )

SES_ALLdivprai.func <- melt(SES_ALLdivprai.func)
colnames(SES_ALLdivprai.func) <- c("diversity","ses")



####Graphs-----
#Now I should be able to plot them in a nice graph
png("Figures/paneltaxandfunct.option2.png", width = 2500, height = 1500, pointsize = 20)
par(mfrow=c(3,4)) # indicates 3 rows, 4 columns
par(mar = c(5,11,4,2)) # sets the margins around the figure

boxplot(ses ~ diversity, data = SES_ALLdivcontrol.tax, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Vacant Lot Taxonomic Diversity", ylim=c(-27,22.5))
stripchart(ses ~ diversity, data = SES_ALLdivcontrol.tax, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(18, 3.2, "A", pos = 4, font = 2, cex = 2.6)

boxplot(ses ~ diversity, data = SES_ALLdivcontrol.func, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Vacant Lot Functional Diversity", ylim=c(-10,10))
stripchart(ses ~ diversity, data = SES_ALLdivcontrol.func, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(8, 4, "B", pos = 4, font = 2, cex = 2.6)

boxplot(ses ~ diversity, data = SES_ALLdivfarm.tax, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Farm Taxonomic Diversity", ylim=c(-27,22.5))
stripchart(ses ~ diversity, data = SES_ALLdivfarm.tax, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(18, 3.2, "C", pos = 4, font = 2, cex = 2.6)

boxplot(ses ~ diversity, data = SES_ALLdivfarm.func, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Farm Functional Diversity", ylim=c(-10,10))
stripchart(ses ~ diversity, data = SES_ALLdivfarm.func, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(8, 4, "D", pos = 4, font = 2, cex = 2.6)




boxplot(ses ~ diversity, data = SES_ALLdivt1.tax, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Vacant Lot Taxonomic Diversity", ylim=c(-27,22.5))
stripchart(ses ~ diversity, data = SES_ALLdivt1.tax, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(18, 3.2, "E", pos = 4, font = 2, cex = 2.6)

boxplot(ses ~ diversity, data = SES_ALLdivt1.func, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Vacant Lot Functional Diversity", ylim=c(-10,10))
stripchart(ses ~ diversity, data = SES_ALLdivt1.func, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(8, 4, "F", pos = 4, font = 2, cex = 2.6)


boxplot(ses ~ diversity, data = SES_ALLdivprai.tax, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Prairie Taxonomic Diversity", ylim=c(-27,22.5))
stripchart(ses ~ diversity, data = SES_ALLdivprai.tax, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(18, 3.2, "G", pos = 4, font = 2, cex = 2.6)

boxplot(ses ~ diversity, data = SES_ALLdivprai.func, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Prairie Functional Diversity", ylim=c(-10, 10))
stripchart(ses ~ diversity, data = SES_ALLdivprai.func, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(8, 4, "H", pos = 4, font = 2, cex = 2.6)


boxplot(ses ~ diversity, data = SES_ALLdivvl.tax, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "2019 Vacant Lot Taxonomic Diversity", ylim=c(-27,22.5))
stripchart(ses ~ diversity, data = SES_ALLdivvl.tax, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(18, 3.2, "I", pos = 4, font = 2, cex = 2.6)

boxplot(ses ~ diversity, data = SES_ALLdivvl.func, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "2019 Vacant Lot Functional Diversity", ylim=c(-10, 10))
stripchart(ses ~ diversity, data = SES_ALLdivvl.func, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(8, 4, "J", pos = 4, font = 2, cex = 2.6)


boxplot(ses ~ diversity, data = SES_ALLdivup.tax, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "2019 Pocket Prairie Taxonomic Diversity", ylim=c(-27,22.5))
stripchart(ses ~ diversity, data = SES_ALLdivup.tax, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(18, 3.2, "K", pos = 4, font = 2, cex = 2.6)


boxplot(ses ~ diversity, data = SES_ALLdivup.func, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "2019 Pocket Prairie Functional Diversity", ylim=c(-10, 10))
stripchart(ses ~ diversity, data = SES_ALLdivup.func, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(8, 4, "L", pos = 4, font = 2, cex = 2.6)

dev.off()


###Making datasets for graphs----
#Okay so to put the tax info from two habitats on one graph and the functional info from two habitats on one graph. Should be just selective copying and pasting from option 2

SES_ALLdiv.mp.tax <- as.data.frame(rbind(UP2019$SES_bsim, UP2019$SES_bsne, UP2019$SES_bsor, VL2019$SES_bsim, VL2019$SES_bsne, VL2019$SES_bsor))
str(SES_ALLdiv.mp.tax)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdiv.mp.tax <- data.frame(t(SES_ALLdiv.mp.tax))
str(SES_ALLdiv.mp.tax)
colnames(SES_ALLdiv.mp.tax) <- c("PP Turnover", "PP Nestedness", " PP Beta-Diversity", "VL Turnover", "VL Nestedness", " VL Beta-Diversity" )


SES_ALLdiv.mp.tax <- melt(SES_ALLdiv.mp.tax)
colnames(SES_ALLdiv.mp.tax) <- c("Taxonomic_diversity","ses")


SES_ALLdiv.mp.func <- as.data.frame(rbind(UP2019$SES_fbsim, UP2019$SES_fbsne, UP2019$SES_fbsor, UP2019$SES_falpha, VL2019$SES_fbsim, VL2019$SES_fbsne, VL2019$SES_fbsor, VL2019$SES_falpha))
str(SES_ALLdiv.mp.func)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdiv.mp.func  <- data.frame(t(SES_ALLdiv.mp.func ))
str(SES_ALLdiv.mp.func )
colnames(SES_ALLdiv.mp.func ) <- c( "PP Turnover", "PP Nestedness",  " PP Beta-Diversity", "PP Alpha", "VL Turnover", "VL Nestedness",  " VL Beta-Diversity", "VL Alpha" )


SES_ALLdiv.mp.func  <- melt(SES_ALLdiv.mp.func )
colnames(SES_ALLdiv.mp.func) <- c("Functional_diversity","ses")


####sivakoff datasets----
SES_ALLdiv.fs.tax <- as.data.frame(rbind(farm$SES_bsim, farm$SES_bsne, farm$SES_bsor, control$SES_bsim, control$SES_bsne, control$SES_bsor))
str(SES_ALLdiv.fs.tax)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdiv.fs.tax <- data.frame(t(SES_ALLdiv.fs.tax))
str(SES_ALLdiv.fs.tax)
colnames(SES_ALLdiv.fs.tax) <- c("UF Turnover", "UF Nestedness", " UF Beta-Diversity", "VL Turnover", "VL Nestedness", " VL Beta-Diversity" )


SES_ALLdiv.fs.tax <- melt(SES_ALLdiv.fs.tax)
colnames(SES_ALLdiv.fs.tax) <- c("Taxonomic_diversity","ses")


SES_ALLdiv.fs.func <- as.data.frame(rbind(farm$SES_fbsim, farm$SES_fbsne, farm$SES_fbsor, farm$SES_falpha, control$SES_fbsim, control$SES_fbsne, control$SES_fbsor, control$SES_falpha))
str(SES_ALLdiv.fs.func)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdiv.fs.func <- data.frame(t(SES_ALLdiv.fs.func))
str(SES_ALLdiv.fs.func)
colnames(SES_ALLdiv.fs.func) <- c( "UF Turnover", "UF Nestedness",  " UF Beta-Diversity", "UF Alpha", "VL Turnover", "VL Nestedness",  " VL Beta-Diversity", "VL Alpha"  )


SES_ALLdiv.fs.func <- melt(SES_ALLdiv.fs.func)
colnames(SES_ALLdiv.fs.func) <- c("Functional_diversity","ses")

####turo datasets----
SES_ALLdiv.kt.tax <- as.data.frame(rbind(Prairie$SES_bsim, Prairie$SES_bsne, Prairie$SES_bsor, T1$SES_bsim, T1$SES_bsne, T1$SES_bsor))
str(SES_ALLdiv.kt.tax)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdiv.kt.tax <- data.frame(t(SES_ALLdiv.kt.tax))
str(SES_ALLdiv.kt.tax)
colnames(SES_ALLdiv.kt.tax) <- c("PP Turnover", "PP Nestedness", " PP Beta-Diversity", "VL Turnover", "VL Nestedness", " VL Beta-Diversity" )


SES_ALLdiv.kt.tax <- melt(SES_ALLdiv.kt.tax)
colnames(SES_ALLdiv.kt.tax) <- c("Taxonomic_diversity","ses")


SES_ALLdiv.kt.func <- as.data.frame(rbind( Prairie$SES_fbsim, Prairie$SES_fbsne, Prairie$SES_fbsor, Prairie$SES_falpha, T1$SES_fbsim, T1$SES_fbsne, T1$SES_fbsor,  T1$SES_falpha))
str(SES_ALLdiv.kt.func)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdiv.kt.func <- data.frame(t(SES_ALLdiv.kt.func))
str(SES_ALLdiv.kt.func)
colnames(SES_ALLdiv.kt.func) <- c( "PP Turnover", "PP Nestedness",  " PP Beta-Diversity", "PP Alpha", "VL Turnover", "VL Nestedness",  " VL Beta-Diversity",  "VL Alpha" )


SES_ALLdiv.kt.func <- melt(SES_ALLdiv.kt.func)
colnames(SES_ALLdiv.kt.func) <- c("Functional_diversity","ses")


####Graphs using for the paper-----
#Going to make three graphs and put them together nicely. 


png("Figures/figure6 panel allfs taxandfunct.png", width = 2500, height = 800, pointsize = 20)
par(mfrow=c(1,2)) # indicates 1 rows, 2 columns
par(mar = c(2,11,2,2)) # sets the margins around the figure (I made them small since I'm going to be editing all the things later)
#Note that to find the colors I used this website (https://contrastchecker.com/) and aimed for a ratio of 2.3 color contrast between the points and boxplot. 
#This website can also be used to find colors with decimal points that will work if you type reb(rdec,gdec,bdec) (https://rgbcolorpicker.com/0-1)
sivakoffpoints.taxcolors<-ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)=="UF Turnover", "#a62b42", 
                       ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)=="UF Nestedness", "#a62b42",
                              ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)==" UF Beta-Diversity", "#a62b42",
                                     ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)=="VL Turnover", "#d57171",
                                            ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)==" VL Beta-Diversity", "#d57171",
                                                   ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)=="VL Nestedness", "#d57171",
                                                          "grey90"))))))
sivakoff.taxcolors<-ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)=="UF Turnover", "#ff5980", 
                             ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)=="UF Nestedness", "#ff5980",
                                    ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)==" UF Beta-Diversity", "#ff5980",
                                           ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)=="VL Turnover", "#ffcccc",
                                                  ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)==" VL Beta-Diversity","#ffcccc",
                                                         ifelse(levels(SES_ALLdiv.fs.tax$Taxonomic_diversity)=="VL Nestedness", "#ffcccc",
                                                                "grey90"))))))
sivakoffpoints.functcolors<-ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="UF Turnover","#a62b42", 
                         ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="UF Nestedness", "#a62b42",
                                ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)==" UF Beta-Diversity", "#a62b42",
                                       ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="UF Alpha", "#a62b42",
                                              ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="VL Turnover", "#d57171",
                                                     ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)==" VL Beta-Diversity", "#d57171",
                                                            ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="VL Alpha", "#d57171",
                                                                   ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="VL Nestedness", "#d57171",
                                                                          "grey90"))))))))
sivakoff.functcolors<-ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="UF Turnover", "#ff5980", 
                               ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="UF Nestedness", "#ff5980",
                                      ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)==" UF Beta-Diversity", "#ff5980",
                                             ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="UF Alpha", "#ff5980",
                                                    ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="VL Turnover", "#ffcccc",
                                                           ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)==" VL Beta-Diversity", "#ffcccc",
                                                                  ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="VL Alpha", "#ffcccc",
                                                                         ifelse(levels(SES_ALLdiv.fs.func$Functional_diversity)=="VL Nestedness", "#ffcccc",
                                                                                "grey90"))))))))




boxplot(ses ~ Taxonomic_diversity, data = SES_ALLdiv.fs.tax, col = sivakoff.taxcolors,
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-32,28), cex.lab = 2, cex.axis=1.45)
stripchart(ses ~ Taxonomic_diversity, data = SES_ALLdiv.fs.tax, col = sivakoffpoints.taxcolors,
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(19.5, 6, "A", pos = 4, font = 2, cex = 2.6)

boxplot(ses ~ Functional_diversity, data = SES_ALLdiv.fs.func, col = sivakoff.functcolors,
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-6,7), cex.lab = 2, cex.axis=1.45)
stripchart(ses ~ Functional_diversity, data = SES_ALLdiv.fs.func, col= sivakoffpoints.functcolors,
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(5.5, 8, "B", pos = 4, font = 2, cex = 2.6)

dev.off()

png("Figures/figure6 panel allkt taxandfunct.png", width = 2500, height = 800, pointsize = 20)
par(mfrow=c(1,2)) # indicates 1 rows, 2 columns
par(mar = c(2,11,2,2)) # sets the margins around the figure (I made them small since I'm going to be editing all the things later)


turopoints.taxcolors<-ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)=="PP Turnover", "#000000", 
                       ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)=="PP Nestedness", "#000000",
                              ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)==" PP Beta-Diversity", "#000000",
                                     ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)=="VL Turnover", "#65A6BF",
                                            ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)==" VL Beta-Diversity", "#65A6BF",
                                                   ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)=="VL Nestedness", "#65A6BF",
                                                          "grey90"))))))
turo.taxcolors<-ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)=="PP Turnover", "#484848", 
                       ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)=="PP Nestedness", "#484848",
                              ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)==" PP Beta-Diversity", "#484848",
                                     ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)=="VL Turnover", "#d1f2ff",
                                            ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)==" VL Beta-Diversity", "#d1f2ff",
                                                   ifelse(levels(SES_ALLdiv.kt.tax$Taxonomic_diversity)=="VL Nestedness", "#d1f2ff",
                                                          "grey90"))))))
turopoints.functcolors<-ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="PP Turnover", "#000000" ,
                         ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="PP Nestedness","#000000",
                                ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)==" PP Beta-Diversity", "#000000",
                                       ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="PP Alpha", "#000000",
                                              ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="VL Turnover", "#65A6BF",
                                                     ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)==" VL Beta-Diversity", "#65A6BF",
                                                            ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="VL Alpha","#65A6BF",
                                                                   ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="VL Nestedness", "#65A6BF",
                                                                          "grey90")))))))) 
turo.functcolors<-ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="PP Turnover", "#484848", 
                         ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="PP Nestedness", "#484848",
                                ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)==" PP Beta-Diversity", "#484848",
                                       ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="PP Alpha", "#484848",
                                              ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="VL Turnover", "#d1f2ff",
                                                     ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)==" VL Beta-Diversity", "#d1f2ff",
                                                            ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="VL Alpha", "#d1f2ff",
                                                                   ifelse(levels(SES_ALLdiv.kt.func$Functional_diversity)=="VL Nestedness", "#d1f2ff",
                                                                          "grey90"))))))))


boxplot(ses ~ Taxonomic_diversity, data = SES_ALLdiv.kt.tax, col = turo.taxcolors,
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, ylim=c(-32,28), cex.lab = 2, cex.axis=1.45)
stripchart(ses ~ Taxonomic_diversity, data = SES_ALLdiv.kt.tax, col = turopoints.taxcolors,
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(19.5, 6, "C", pos = 4, font = 2, cex = 2.6)

boxplot(ses ~ Functional_diversity, data = SES_ALLdiv.kt.func, col = turo.functcolors,
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-6,7), cex.lab = 2, cex.axis=1.45)
stripchart(ses ~ Functional_diversity, data = SES_ALLdiv.kt.func, col = turopoints.functcolors,
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(5.5, 8, "D", pos = 4, font = 2, cex = 2.6)
dev.off()

png("Figures/figure6 panel allmp taxandfunct.png", width = 2500, height = 800, pointsize = 20)
par(mfrow=c(1,2)) # indicates 1 rows, 2 columns
par(mar = c(2,11,2,2)) # sets the margins around the figure (I made them small since I'm going to be editing all the things later)


phampoints.taxcolors<-ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)=="PP Turnover", "#5e0a7c", 
                       ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)=="PP Nestedness", "#5e0a7c",
                              ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)==" PP Beta-Diversity", "#5e0a7c",
                                     ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)=="VL Turnover", "#985dc7",
                                     ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)==" VL Beta-Diversity", "#985dc7",
                                            ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)=="VL Nestedness", "#985dc7",
                                                   "grey90"))))))
pham.taxcolors<-ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)=="PP Turnover","#a340cc", 
                       ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)=="PP Nestedness", "#a340cc",
                              ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)==" PP Beta-Diversity", "#a340cc",
                                     ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)=="VL Turnover", "#d9a6ff",
                                            ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)==" VL Beta-Diversity", "#d9a6ff",
                                                   ifelse(levels(SES_ALLdiv.mp.tax$Taxonomic_diversity)=="VL Nestedness", "#d9a6ff",
                                                          "grey90"))))))
phampoints.functcolors<-ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="PP Turnover", "#5e0a7c", 
                       ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="PP Nestedness", "#5e0a7c",
                              ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)==" PP Beta-Diversity", "#5e0a7c",
                                     ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="PP Alpha", "#5e0a7c",
                                     ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="VL Turnover", "#985dc7",
                                            ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)==" VL Beta-Diversity", "#985dc7",
                                                   ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="VL Alpha", "#985dc7",
                                                   ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="VL Nestedness", "#985dc7",
                                                          "grey90"))))))))

pham.functcolors<-ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="PP Turnover", "#a340cc", 
                         ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="PP Nestedness", "#a340cc",
                                ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)==" PP Beta-Diversity", "#a340cc",
                                       ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="PP Alpha", "#a340cc",
                                              ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="VL Turnover", "#d9a6ff",
                                                     ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)==" VL Beta-Diversity", "#d9a6ff",
                                                            ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="VL Alpha", "#d9a6ff",
                                                                   ifelse(levels(SES_ALLdiv.mp.func$Functional_diversity)=="VL Nestedness", "#d9a6ff",
                                                                          "grey90"))))))))

boxplot(ses ~ Taxonomic_diversity, data = SES_ALLdiv.mp.tax, col = pham.taxcolors,
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-32,28), cex.lab = 2, cex.axis=1.45)
stripchart(ses ~ Taxonomic_diversity, data = SES_ALLdiv.mp.tax, col = phampoints.taxcolors,
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(19.5, 6, "E", pos = 4, font = 2, cex = 2.6)

boxplot(ses ~ Functional_diversity, data = SES_ALLdiv.mp.func, col = pham.functcolors,
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-6,7), cex.lab = 2, cex.axis=1.45)
stripchart(ses ~ Functional_diversity, data = SES_ALLdiv.mp.func, col = phampoints.functcolors,
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(5.5, 8, "F", pos = 4, font = 2, cex = 2.6)



dev.off()






##Figure of Turo sig figs----
library(ggplot2)
library(cowplot)
treatmentcleanup.mp <- c(VL="VL", UP="PP")
SES.mp$treatmentspelledout <- as.character(treatmentcleanup.mp[SES.mp$trmt])
treatmentcleanup.kt <- c(T1="VL", Prairie="PP")
SES.kt$treatmentspelledout <- as.character(treatmentcleanup.kt[SES.kt$trmt])



png("Figures/fig7 panel turophamsigs.png", width = 1500, height = 1000, pointsize = 20)
windows()
par(mfrow=c(1,2)) # indicates one row, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure

# Functional Alpha
boxplot(SES_falpha ~ treatmentspelledout, data = SES.kt, col = c("#C1BDFF", "#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_falpha ~ trmt, data = SES.kt, col = c("#5e0a7c","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.87, 2.42, "A", pos = 4, font = 2, cex = 2)

# Soil Nesting.kt
boxplot(SES_nest_1 ~ treatmentspelledout, data = SES.kt,  col = c("#C1BDFF", "#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_nest_1 ~ trmt, data = SES.kt, col = c("#5e0a7c","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.87, 2.42, "B", pos = 4, font = 2, cex = 2)

# Parasitic
#boxplot(SES_soc_4 ~ treatmentspelledout, data = SES.mp, col = c("#C1BDFF", "#52A43B"),
#        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
#        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
#        horizontal = TRUE, las = 1, range = 0)
#stripchart(SES_soc_4 ~ trmt, data = SES.mp, col = c("#5e0a7c","#11290A"),
#           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
#abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.87, 2.42, "C", pos = 4, font = 2, cex = 2)

# Soil Nesting mp
#boxplot(SES_nest_1 ~ treatmentspelledout, data = SES.mp, col = c("#C1BDFF", "#52A43B"),
#        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
#        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
#        horizontal = TRUE, las = 1, range = 0 )
#stripchart(SES_nest_1 ~ trmt, data = SES.mp, col = c("#5e0a7c","#11290A"),
 #          pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
#abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.87, 2.42, "D", pos = 4, font = 2, cex = 2)



dev.off()


##Figure of sivakoff sig figs----




treatmentcleanup.fs <- c(Farm="UF", Control="VL")
SES.fs$treatmentspelledout <- as.character(treatmentcleanup.fs[SES.fs$trmt])


png("Figures/fig 6 panel sivakoffsigs.png", width = 1500, height = 1000, pointsize = 20)

par(mfrow=c(2,2)) # indicates one row, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure



# Soil Nesting.kt
boxplot(SES_nest_1 ~ treatmentspelledout, data = SES.fs,  col = c("#6D3737","#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_nest_1 ~ treatmentspelledout, data = SES.fs, col = c("burlywood2","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.87, 2.42, "B", pos = 4, font = 2, cex = 2)

# Cavity nesting
boxplot(SES_nest_2 ~ treatmentspelledout, data = SES.fs, col = c("#6D3737","#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45,
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_nest_2 ~ treatmentspelledout, data = SES.fs, col = c("burlywood2","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(2.38, 2.41, "B", pos = 4, font = 2, cex = 2)

# native
boxplot(SES_ori_0 ~ treatmentspelledout, data = SES.fs, col = c("#6D3737","#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45,
        horizontal = TRUE, las = 1, range = 0 )
stripchart(SES_ori_0 ~ treatmentspelledout, data = SES.fs, col = c("burlywood2","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(2.38, 2.41, "C", pos = 4, font = 2, cex = 2)


# non-native
boxplot(SES_ori_1 ~ treatmentspelledout, data = SES.fs, col =c("#6D3737","#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45,
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_ori_1 ~ treatmentspelledout, data = SES.fs, col = c("burlywood2","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(2.38, 2.41, "D", pos = 4, font = 2, cex = 2)


dev.off()

#My favorite website for picking out colors: https://www.whocanuse.com/
#Figure of pham traits----
png("Figures/fig panel phamsigs.png", width = 1500, height = 1000, pointsize = 20)

par(mfrow=c(2,2)) # indicates two rows, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure

# Functional Alpha
boxplot(SES_falpha ~ treatmentspelledout, data = SES.mp, col = c("#C1BDFF", "#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_falpha ~ trmt, data = SES.mp, col = c("#5e0a7c","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.87, 2.42, "A", pos = 4, font = 2, cex = 2)

# Soil Nesting.mp
boxplot(SES_nest_1 ~ treatmentspelledout, data = SES.mp,  col = c("#C1BDFF", "#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_nest_1 ~ trmt, data = SES.mp, col = c("#5e0a7c","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.87, 2.42, "B", pos = 4, font = 2, cex = 2)

# Parasitic
#boxplot(SES_soc_4 ~ treatmentspelledout, data = SES.mp, col = c("#C1BDFF", "#52A43B"),
#        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
#        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
#        horizontal = TRUE, las = 1, range = 0)
#stripchart(SES_soc_4 ~ trmt, data = SES.mp, col = c("#5e0a7c","#11290A"),
#           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
#abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.87, 2.42, "C", pos = 4, font = 2, cex = 2)

# Soil Nesting mp
#boxplot(SES_nest_1 ~ treatmentspelledout, data = SES.mp, col = c("#C1BDFF", "#52A43B"),
#        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
#        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
#        horizontal = TRUE, las = 1, range = 0 )
#stripchart(SES_nest_1 ~ trmt, data = SES.mp, col = c("#5e0a7c","#11290A"),
#          pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
#abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(1.87, 2.42, "D", pos = 4, font = 2, cex = 2)



dev.off()



#New.Graphs. Pham and TuroFunctional diversities----
##Turo funct dataset----
SES_ALLdiv.kt.func.forcomp <- as.data.frame(rbind( Prairie$SES_fbsim, T1$SES_fbsim, Prairie$SES_fbsne,T1$SES_fbsne, Prairie$SES_fbsor, T1$SES_fbsor, Prairie$SES_falpha, T1$SES_falpha))
str(SES_ALLdiv.kt.func.forcomp)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdiv.kt.func.forcomp <- data.frame(t(SES_ALLdiv.kt.func.forcomp))
str(SES_ALLdiv.kt.func.forcomp)
colnames(SES_ALLdiv.kt.func.forcomp) <- c( "PP Turnover", "VL Turnover", "PP Nestedness", "VL Nestedness",  " PP Beta-Diversity", " VL Beta-Diversity",  "PP Alpha",  "VL Alpha" )


SES_ALLdiv.kt.func.forcomp <- melt(SES_ALLdiv.kt.func.forcomp)
colnames(SES_ALLdiv.kt.func.forcomp) <- c("Functional_diversity","ses")

##Pham funct dataset----

SES_ALLdiv.mp.func.forcomp <- as.data.frame(rbind(UP2019$SES_fbsim, VL2019$SES_fbsim, UP2019$SES_fbsne, VL2019$SES_fbsne, UP2019$SES_fbsor, VL2019$SES_fbsor, UP2019$SES_falpha,  VL2019$SES_falpha))
str(SES_ALLdiv.mp.func.forcomp)
#Not sure where that last column comes from
#It bound them in the wrong direction. The t function should switch row and columns. 
SES_ALLdiv.mp.func.forcomp  <- data.frame(t(SES_ALLdiv.mp.func.forcomp ))
str(SES_ALLdiv.mp.func.forcomp )
colnames(SES_ALLdiv.mp.func.forcomp ) <- c( "PP Turnover", "VL Turnover", "PP Nestedness", "VL Nestedness",  " PP Beta-Diversity", " VL Beta-Diversity",  "PP Alpha",  "VL Alpha"  )


SES_ALLdiv.mp.func.forcomp  <- melt(SES_ALLdiv.mp.func.forcomp )
colnames(SES_ALLdiv.mp.func.forcomp) <- c("Functional_diversity","ses")

##colors----
turo.functcolors.forcomp<-ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="PP Turnover", "#C1BDFF" ,
                               ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="PP Nestedness","#C1BDFF",
                                      ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)==" PP Beta-Diversity", "#C1BDFF",
                                             ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="PP Alpha", "#C1BDFF",
                                                    ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="VL Turnover", "#52A43B",
                                                           ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)==" VL Beta-Diversity", "#52A43B",
                                                                  ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="VL Alpha","#52A43B",
                                                                         ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="VL Nestedness", "#52A43B",
                                                                                "grey90")))))))) 
turopoints.functcolors.forcomp <-ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="PP Turnover", "#5e0a7c", 
                         ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="PP Nestedness", "#5e0a7c",
                                ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)==" PP Beta-Diversity", "#5e0a7c",
                                       ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="PP Alpha", "#5e0a7c",
                                              ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="VL Turnover", "#11290A",
                                                     ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)==" VL Beta-Diversity", "#11290A",
                                                            ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="VL Alpha", "#11290A",
                                                                   ifelse(levels(SES_ALLdiv.kt.func.forcomp$Functional_diversity)=="VL Nestedness", "#11290A",
                                                                          "grey90"))))))))

pham.functcolors.forcomp<-ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="PP Turnover", "#C1BDFF", 
                               ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="PP Nestedness", "#C1BDFF",
                                      ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)==" PP Beta-Diversity", "#C1BDFF",
                                             ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="PP Alpha", "#C1BDFF",
                                                    ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="VL Turnover", "#52A43B",
                                                           ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)==" VL Beta-Diversity", "#52A43B",
                                                                  ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="VL Alpha", "#52A43B",
                                                                         ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="VL Nestedness", "#52A43B",
                                                                                "grey90"))))))))

phampoints.functcolors.forcomp <-ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="PP Turnover", "#5e0a7c", 
                         ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="PP Nestedness", "#5e0a7c",
                                ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)==" PP Beta-Diversity", "#5e0a7c",
                                       ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="PP Alpha", "#5e0a7c",
                                              ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="VL Turnover", "#11290A",
                                                     ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)==" VL Beta-Diversity", "#11290A",
                                                            ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="VL Alpha", "#11290A",
                                                                   ifelse(levels(SES_ALLdiv.mp.func.forcomp$Functional_diversity)=="VL Nestedness", "#11290A",
                                                                          "grey90"))))))))


##graphs----
png("Figures/ panel ktvsmp functdiv.png", width = 2500, height = 800, pointsize = 20)
par(mfrow=c(1,2)) # indicates 1 rows, 2 columns
par(mar = c(2,11,2,2)) # sets the margins around the figure (I made them small since I'm going to be editing all the things later)

boxplot(ses ~ Functional_diversity, data = SES_ALLdiv.kt.func.forcomp, col = turo.functcolors.forcomp,
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-6,7), cex.lab = 2, cex.axis=1.45)
stripchart(ses ~ Functional_diversity, data = SES_ALLdiv.kt.func.forcomp, col = turopoints.functcolors.forcomp,
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
#text(5.5, 8, "D", pos = 4, font = 2, cex = 2.6)

boxplot(ses ~ Functional_diversity, data = SES_ALLdiv.mp.func.forcomp, col = pham.functcolors.forcomp,
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0,  ylim=c(-6,7), cex.lab = 2, cex.axis=1.45)
stripchart(ses ~ Functional_diversity, data = SES_ALLdiv.mp.func.forcomp, col = phampoints.functcolors.forcomp,
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

dev.off()


png("Figures/fig7  turo and pham soil.png", width = 1500, height = 1000, pointsize = 20)
par(mfrow=c(2,2)) # indicates one row, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure
# Soil Nesting.kt
boxplot(SES_nest_1 ~ treatmentspelledout, data = SES.kt,  col = c("#C1BDFF", "#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_nest_1 ~ trmt, data = SES.kt, col = c("#5e0a7c","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

# Soil Nesting.mp
boxplot(SES_nest_1 ~ treatmentspelledout, data = SES.mp,  col = c("#C1BDFF", "#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_nest_1 ~ trmt, data = SES.mp, col = c("#5e0a7c","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

dev.off()

png("Figures/fig7b  turo and pham pithy stem.png", width = 1500, height = 1000, pointsize = 20)
par(mfrow=c(2,2)) # indicates one row, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure
# pithy stem Nesting.kt
boxplot(SES_nest_4 ~ treatmentspelledout, data = SES.kt,  col = c("#C1BDFF", "#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_nest_4 ~ trmt, data = SES.kt, col = c("#5e0a7c","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

# pithy stem Nesting.mp
boxplot(SES_nest_4 ~ treatmentspelledout, data = SES.mp,  col = c("#C1BDFF", "#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_nest_4 ~ trmt, data = SES.mp, col = c("#5e0a7c","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

dev.off()

png("Figures/fig7c  turo and pham solitary.png", width = 1500, height = 1000, pointsize = 20)
par(mfrow=c(2,2)) # indicates one row, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure
# solitary.kt
boxplot(SES_soc_2 ~ treatmentspelledout, data = SES.kt,  col = c("#C1BDFF", "#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_soc_2 ~ trmt, data = SES.kt, col = c("#5e0a7c","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

# solitary.mp
boxplot(SES_soc_2 ~ treatmentspelledout, data = SES.mp,  col = c("#C1BDFF", "#52A43B"),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-3.6,3), cex.lab = 2, cex.axis = 1.45, 
        horizontal = TRUE, las = 1, range = 0)
stripchart(SES_soc_2 ~ trmt, data = SES.mp, col = c("#5e0a7c","#11290A"),
           pch = 19, cex = 1.2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

dev.off()

