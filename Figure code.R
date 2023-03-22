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
