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
stripchart(SES_ori_0 ~ trmt, data = fc, col = viridis(3),
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
stripchart(SES_ori_0 ~ trmt, data = fc, col = viridis(3),
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


library(ggplot2)





#Length and Origin
SES_lengthandorigin.FC<- data.frame(
  trmt=rep(fc$trmt,3),
  variable=rep(c("Body Length", "Native", "Non-Native"), each = nrow(fc)),
  value = c(fc$SES_bl, fc$SES_ori_0, fc$SES_ori_1),
  ordered=TRUE)


png("Figure 5A all FS local sp traits.png", width = 1500, height = 1000, pointsize = 20)

ggplot(SES_lengthandorigin.FC, aes(x = factor(variable,levels = c("Body Length", "Native", "Non-Native"), ordered=TRUE), y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             shape = 19, size = 5, show.legend = FALSE) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=25)) +
  theme(axis.text.y = element_text(size=39))+
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip() +guides(fill = FALSE)


dev.off()

#Nesting traits

SES_nesting.FC<- data.frame(
  trmt=rep(fc$trmt,5),
  variable=rep(c("Soil","Cavity", "Colony", "Pithy Stems", "Wood"), each = nrow(fc)),
  value = c(fc$SES_nest_1, fc$SES_nest_2, fc$SES_nest_3, fc$SES_nest_4, fc$SES_nest_5)
)

png("Figure 5B all FS local sp traits.png", width = 1500, height = 1000, pointsize = 20)

ggplot(SES_nesting.FC, aes(x = variable, y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             shape = 19, size = 5, show.legend = FALSE) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size=25)) +
  theme(axis.text.y = element_text(size=39))+
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip()+ guides(fill = FALSE)
dev.off()

#Lecty
SES_lecty.FC<- data.frame(
  trmt=rep(fc$trmt,3),
  variable=rep(c("Kleptoparasitic","Generalist", "Specialist"), each = nrow(fc)),
  value = c(fc$SES_lec_0, fc$SES_lec_1, fc$SES_lec_2)
)

png("Figure 5C all FS local sp traits.png", width = 1500, height = 1000, pointsize = 20)

 ggplot(SES_lecty.FC, aes(x = variable, y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             shape = 19, size = 5, show.legend = FALSE) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=25)) +
   theme(axis.text.y = element_text(size=39))+
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip()+ guides(fill = FALSE)
dev.off()

#Sociality
SES_sociality.FC<- data.frame(
  trmt=rep(fc$trmt,4),
  variable=rep(c("Solitary","Subsocial", "Eusocial", "Parasitic"), each = nrow(fc)),
  value = c(fc$SES_soc_1, fc$SES_soc_2, fc$SES_soc_3, fc$SES_soc_4)
)

png("Figure 5D all FS local sp traits.png", width = 1500, height = 1000, pointsize = 20)

ggplot(SES_sociality.FC, aes(x = variable, y = value, fill = trmt)) +
  geom_boxplot(position = position_dodge(width = .75), alpha = 0.6, coef = 0, width = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = .75), 
             shape = 19, size = 5, show.legend = FALSE) +
  scale_fill_viridis_d(option = "C", end = 0.3, direction = -1, alpha = 0.6) +
  ylab("Standardized Effect Sizes (SES)") + 
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=25)) +
  theme(axis.text.y = element_text(size=39))+
  ylim(c(-8, 8)) +
  geom_hline(yintercept = 0, col = "black", lwd = 2, linetype = "dashed") +
  coord_flip() +guides(fill = FALSE)



dev.off()




