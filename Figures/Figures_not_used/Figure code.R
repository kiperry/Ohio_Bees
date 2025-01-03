##################################################################################################
# from line 452 in ComparingLocal file

# Make figures!

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
png("Figures/Figure 1.png", width = 1500, height = 1000, pointsize = 20)

par(mfrow=c(2,2)) # indicates two rows, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure

# hive nesting
boxplot(SES_nest_3 ~ trmt, data = fc, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-2,3), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Hive Nesting")
stripchart(SES_nest_3 ~ trmt, data = fc, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

# native
boxplot(SES_ori_0 ~ trmt, data = fc, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-2,3), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Native")
stripchart(SES_ori_0 ~ trmt, data = fc, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

# specialists
boxplot(SES_lec_2 ~ trmt, data = fc, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "", 
        ylim = c(-2,3), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Specialists")
stripchart(SES_lec_2 ~ trmt, data = fc, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

# non-native
boxplot(SES_ori_1 ~ trmt, data = fc, col = viridis(3, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        ylim = c(-2,3), cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.5,
        horizontal = TRUE, las = 1, range = 0, main = "Non-Native")
stripchart(SES_ori_1 ~ trmt, data = fc, col = viridis(3),
           pch = 19, cex = 2, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)


dev.off()

##################################################################################################

#Okay Carlee's attempt to make figures
#Figure to replace table 2. Want it to look like 2 panels one of the vacant lots one of the Urban farms
#In each panel is the taxonomic beta, functional alpha functional beta and all components listed in text. 
#order: functional alpha, taxonomic beta, (turnover , nestedness), funct. beta (T,N)


##Loading programs needed
if (!suppressWarnings(require(viridis))) install.packages("viridis")
citation("viridis")

if (!suppressWarnings(require(reshape2))) install.packages("reshape2")
citation("reshape2")


#step one- create a dataframe for all of the tax and funct diversity metrics I want to graph
SES_ALLdivcontrol <- as.data.frame(rbind( control$SES_fbsim, control$SES_fbsne, control$SES_fbsor, control$SES_bsim, control$SES_bsne, control$SES_bsor, control$SES_falpha))
str(SES_ALLdivcontrol)
#Not sure where that last column comes from
#OHHH it bound them in the wrong direction. The t function should switch row and columns. 
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
png("Local Diversity Metrics Figure.png", width = 1500, height = 1000, pointsize = 20)
par(mfrow=c(1,2)) # indicates one row, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure

boxplot(ses ~ diversity, data = SES_ALLdivcontrol, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Vacant Lot Diversity")
stripchart(ses ~ diversity, data = SES_ALLdivcontrol, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

boxplot(ses ~ diversity, data = SES_ALLdivfarm, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Farm Diversity")
stripchart(ses ~ diversity, data = SES_ALLdivfarm, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
dev.off()

#Okay I don't know how to fix the Words, but the graphs exist!!!! Saving and uploadign to git


#Next we want to make a figure that has the four panels, of lenght and origin, nesting traits, lecty and soicality
#But with both urban farms and vacant lots. 


#Below is one attempt to make that figure. Will determine if a different version is preferable

library(ggplot2)
library(cowplot)




#Length and Origin
SES_lengthandorigin.FC<- data.frame(
  trmt=rep(fc$trmt,3),
  variable=rep(c("Alien", "Native", "Body Length"), each = nrow(fc)),
  value = c(fc$SES_ori_1, fc$SES_ori_0, fc$SES_bl)
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
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip() +guides(fill = FALSE)



#Nesting traits

SES_nesting.FC<- data.frame(
  trmt=rep(fc$trmt,5),
  variable=rep(c("Wood","Pithy Stems", "Colony", "Cavity", "Soil"), each = nrow(fc)),
  value = c(fc$SES_nest_5, fc$SES_nest_4, fc$SES_nest_3, fc$SES_nest_2, fc$SES_nest_1)
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
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip()+ guides(fill = FALSE)


#Lecty
SES_lecty.FC<- data.frame(
  trmt=rep(fc$trmt,3),
  variable=rep(c("Specialist","Generalist", "Kleptoparasitic"), each = nrow(fc)),
  value = c(fc$SES_lec_2, fc$SES_lec_1, fc$SES_lec_0)
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
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip()+ guides(fill = FALSE)


#Sociality
SES_sociality.FC<- data.frame(
  trmt=rep(fc$trmt,4),
  variable=rep(c("Parasitic", "Eusocial","Subsocial", "Solitary"), each = nrow(fc)),
  value = c(fc$SES_soc_4, fc$SES_soc_3, fc$SES_soc_2, fc$SES_soc_1)
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
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip() +guides(fill = FALSE)


fig5 <- plot_grid(fig5a,fig5b,fig5c,fig5d, labels= c('A','B','C','D'))
png("Figure 5 all FS local sp traits.png", width = 1500, height = 1000, pointsize = 20)

fig5

dev.off()

#The below code is used in tandem with the ComparingLocal.R file, which had the code for the plots already written


a6<- effect_plot(FS_bl.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = 'Percentage Greenspace', y.label = 'SES', main.title = 'Body Length')
b6<- effect_plot(FS_bl.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Greenspace Isolation',y.label ='', main.title = 'Body Length')
c6<- effect_plot(FS_ori_0.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = ' Greenspace Isolation',y.label ='', main.title = 'Native Species')
d6<- effect_plot(FS_ori_1.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Greenspace Isolation ', y.label ='SES', main.title = 'Alien Species')
e6<- effect_plot(FS_nest_1.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'Treatment', y.label ='', main.title = 'Soil Nesting')
f6<- effect_plot(FS_nest_2.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'Treatment', y.label ='', main.title = 'Cavity Nesting')
g6<- effect_plot(FS_lec_2.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = ' Percentage Greenspace', y.label ='SES', main.title = 'Specialist')


install.packages("cowplot")
library(cowplot)
fig6 <- plot_grid(a6,b6,c6,d6,e6,f6,g6, labels= c('A','B','C','D','E','F','G'))

png("Figure 6 All functional predictions.png", width = 1500, height = 1000, pointsize = 20)

fig6
dev.off()


#Next making figure of vacant lot and pocket prairie diversity metrics

t1t8 <- rbind(T1, T8)
str(t1t8)
t1t8 <- droplevels(t1t8) #drops farm and control
str(t1t8)

# changes the treatment names for the figure
levels(t1t8$trmt)[levels(t1t8$trmt)=='T1'] <- 'Vacant Lot'
levels(t1t8$trmt)[levels(t1t8$trmt)=='T8'] <- 'Pocket Prairie'
levels(t1t8$trmt)

#step one- create a dataframe for all of the tax and funct diversity metrics I want to graph
SES_ALLdivvl <- as.data.frame(rbind( T1$SES_fbsim, T1$SES_fbsne, T1$SES_fbsor, T1$SES_bsim, T1$SES_bsne, T1$SES_bsor, T1$SES_falpha))
str(SES_ALLdivvl)
SES_ALLdivvl <- data.frame(t(SES_ALLdivvl)) #t function switches columns with rows
str(SES_ALLdivvl)
colnames(SES_ALLdivvl) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity","Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity", "Funct. Alpha" )


SES_ALLdivvl <- melt(SES_ALLdivvl)
colnames(SES_ALLdivvl) <- c("diversity","ses")



SES_ALLdivprai <- as.data.frame(rbind( T8$SES_fbsim, T8$SES_fbsne, T8$SES_fbsor, T8$SES_bsim, T8$SES_bsne, T8$SES_bsor, T8$SES_falpha))
str(SES_ALLdivprai)
SES_ALLdivprai <- data.frame(t(SES_ALLdivprai))
str(SES_ALLdivprai)
colnames(SES_ALLdivprai) <- c( "Funct. Turnover", "Funct. Nestedness",  "Total Funct. Beta-Diversity","Tax. Turnover", "Tax. Nestedness", "Total Tax. Beta-Diversity", "Funct. Alpha")

SES_ALLdivprai <- melt(SES_ALLdivprai)
colnames(SES_ALLdivprai) <- c("diversity","ses")


#Now I should be able to plot it in a nice graph
png("Local Diversity Metrics Figure7.png", width = 1500, height = 1000, pointsize = 20)
par(mfrow=c(1,2)) # indicates one row, two columns
par(mar = c(5,7,4,2)) # sets the margins around the figure

boxplot(ses ~ diversity, data = SES_ALLdivvl, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Vacant Lot Diversity")
stripchart(ses ~ diversity, data = SES_ALLdivvl, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)

boxplot(ses ~ diversity, data = SES_ALLdivprai, col = viridis(6, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, main = "Pocket Prairie Diversity")
stripchart(ses ~ diversity, data = SES_ALLdivprai, col = viridis(6),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
dev.off()



#Making Figure 8

library(ggplot2)





#Length and Origin
SES_lengthandorigin.t1t8<- data.frame(
  trmt=rep(t1t8$trmt,3),
  variable=rep(c("Alien", "Native", "Body Length"), each = nrow(t1t8)),
  value = c(t1t8$SES_ori_1, t1t8$SES_ori_0, t1t8$SES_bl)
)
SES_lengthandorigin.t1t8$variable<-factor(SES_lengthandorigin.t1t8$variable, c("Alien", "Native", "Body Length"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots



fig8a<- ggplot(SES_lengthandorigin.t1t8, aes(x = variable, y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             shape = 19, size = 3, show.legend = TRUE) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=15)) +
  theme(axis.text.y = element_text(size=25))+
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip() 



#Nesting traits

SES_nesting.t1t8<- data.frame(
  trmt=rep(t1t8$trmt,5),
  variable=rep(c("Wood","Pithy Stems", "Colony", "Cavity", "Soil"), each = nrow(t1t8)),
  value = c(t1t8$SES_nest_5, t1t8$SES_nest_4, t1t8$SES_nest_3, t1t8$SES_nest_2, t1t8$SES_nest_1)
)
SES_nesting.t1t8$variable<-factor(SES_nesting.t1t8$variable, c("Wood", "Pithy Stems","Colony", "Cavity", "Soil"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots


fig8b<- ggplot(SES_nesting.t1t8, aes(x = variable, y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             shape = 19, size = 3, show.legend = TRUE) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size=15)) +
  theme(axis.text.y = element_text(size=25))+
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip()


#Lecty
SES_lecty.t1t8<- data.frame(
  trmt=rep(t1t8$trmt,3),
  variable=rep(c("Specialist","Generalist", "Kleptoparasitic"), each = nrow(t1t8)),
  value = c(t1t8$SES_lec_2, t1t8$SES_lec_1, t1t8$SES_lec_0)
)
SES_lecty.t1t8$variable<-factor(SES_lecty.t1t8$variable, c("Specialist", "Generalist","Kleptoparasitic"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots


fig8c<- ggplot(SES_lecty.t1t8, aes(x = variable, y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             shape = 19, size = 3, show.legend = TRUE) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=15)) +
  theme(axis.text.y = element_text(size=25))+
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip()


#Sociality
SES_sociality.t1t8<- data.frame(
  trmt=rep(t1t8$trmt,4),
  variable=rep(c("Parasitic", "Eusocial","Subsocial", "Solitary"), each = nrow(t1t8)),
  value = c(t1t8$SES_soc_4, t1t8$SES_soc_3, t1t8$SES_soc_2, t1t8$SES_soc_1)
) #Note that we have to put our funct traits in the reverse order of how we want them to appear on our graph
SES_sociality.t1t8$variable<-factor(SES_sociality.t1t8$variable, c("Parasitic", "Eusocial","Subsocial", "Solitary"))
#The above code puts the funct. traits as "factors" which is necessary for us to keep the order we want in our box plots

fig8d<- ggplot(SES_sociality.t1t8, aes(x = variable, y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = .75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = .75), 
             shape = 19, size = 3) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") + 
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=15)) +
  theme(axis.text.y = element_text(size=25))+
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip() 


fig8 <- plot_grid(fig8a,fig8b,fig8c,fig8d, labels= c('A','B','C','D'))
png("Figure 8 all FS local sp traits.png", width = 1500, height = 1000, pointsize = 20)

fig8

dev.off()

#figure 9



a9<- effect_plot(KT_falpha.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'Treatment', y.label = 'Standardized Effect Sizes (SES)')
b9<- effect_plot(KT_bsor.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Greenspace Isolation (ENN)', y.label = 'Standardized Effect Sizes (SES)')


install.packages("cowplot")
library(cowplot)
fig9 <- plot_grid(a9,b9, labels= c('A','B'))
png("Figures/Figure 9.png", width = 1500, height = 1000, pointsize = 20)

fig9
dev.off()



#fig 10

a10<- effect_plot(KT_nest_3.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'Treatment', y.label = 'Standardized Effect Sizes (SES)', main.title = 'Colony Nesting')
b10<- effect_plot(KT_nest_5.mod.red, pred = trmt, interval = TRUE, partial.residuals = TRUE, x.label = 'Treatment', y.label = 'Standardized Effect Sizes (SES)', main.title = 'Wood Nesting')
c10<- effect_plot(KT_nest_4.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = 'Percentage Greenspace', y.label = 'Standardized Effect Sizes (SES)', main.title = 'Pithy Stem Nesting')
d10<- effect_plot(KT_nest_4.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Greenspace Isolation', y.label = 'Standardized Effect Sizes (SES)', main.title = 'Pithy Stem Nesting')
e10<- effect_plot(KT_soc_1.mod.red, pred = pland, interval = TRUE, partial.residuals = TRUE, x.label = 'Percentage Greenspace', y.label = 'Standardized Effect Sizes (SES)', main.title = 'Subsocial')
f10<- effect_plot(KT_soc_1.mod.red, pred = enn, interval = TRUE, partial.residuals = TRUE, x.label = 'Greenspace Isolation', y.label = 'Standardized Effect Sizes (SES)', main.title = 'Subsocial')

fig10 <- plot_grid(a10,b10,c10,d10,e10,f10, labels= c('A','B','C','D','E','F'))
png("Figures/Figure 10.png", width = 1500, height = 1000, pointsize = 20)

fig10
dev.off()
