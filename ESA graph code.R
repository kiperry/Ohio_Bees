#Let's try making some figures for the local things, shall we?

#loading the things
library(ggplot2)
library(ggthemes)
library(viridis)

#Make sure the local analysis workspace is open. 

#Okay this is the way to make the graph for specialists
boxplot(SES_lec_2 ~ trmt, data = SES, col = viridis(4, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, cex.lab = 1.5, cex.axis = 1)
stripchart(SES_lec_2 ~ trmt, data = SES, col = viridis(4),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(4.5, 4.2, "C", pos = 4, font = 2, cex = 2.6)

#For wood nesting
boxplot(SES_nest_5 ~ trmt, data = SES, col = viridis(4, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, cex.lab = 1.5, cex.axis = 1)
stripchart(SES_nest_5 ~ trmt, data = SES, col = viridis(4),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(4.5, 4.2, "C", pos = 4, font = 2, cex = 2.6)

#For soil nesting
boxplot(SES_nest_1 ~ trmt, data = SES, col = viridis(4, alpha = 0.6),
        xlab = "Standardized Effect Sizes (SES)", ylab = "",
        horizontal = TRUE, las = 1, range = 0, cex.lab = 1.5, cex.axis = 1)
stripchart(SES_nest_1 ~ trmt, data = SES, col = viridis(4),
           pch = 19, cex = 1, las = 1, add = TRUE, method = "jitter", jitter = 0.2)
abline(v = 0.0, col = "black", lwd = 3, lty=2)
text(4.5, 4.2, "c", pos = 4, font = 2, cex = 2.6)



dev.off()
