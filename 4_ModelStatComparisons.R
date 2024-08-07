library(dplyr)
library(ggplot2)
library(ggtext)
library(ggpubr)

# Read in data
glm2D <- read.csv("~/Dropbox/MARDIGRA/Papers/Manuscript/Appendix Material/GLM2DModelDiagnosticStats.csv")
glm3D <- read.csv("~/Dropbox/MARDIGRA/Papers/Manuscript/Appendix Material/GLM3DModelDiagnosticStats.csv")

ogNames <- names(glm3D)

glmDiff <- glm3D[,-1] - glm2D[,-1]

names(glm2D)[-1] <- paste0(names(glm2D)[-1], "_2D")
names(glm3D)[-1] <- paste0(names(glm3D)[-1], "_3D")
names(glmDiff) <- paste0(names(glmDiff), "_3D-2D")
glmDiff$species <- glm2D$species

# Look at relevant correlations
pearCor <- cor(glm2D[,5:10], glm3D[,5:10], 
               use = "pairwise.complete.obs")
write.csv(pearCor, "~/Dropbox/MARDIGRA/Papers/Manuscript/Appendix Material/CorrelationAmongModelStats.csv")

allGLMStats <- merge(glm2D, glm3D)
allGLMStats <- merge(allGLMStats, glmDiff)

# Statistical tests
wilcox.test(allGLMStats$nPres_2D, allGLMStats$nPres_3D, 
            paired=TRUE, alternative = "less")

wilcox.test(allGLMStats$nAbs_2D, allGLMStats$nAbs_3D, 
            paired=TRUE, alternative = "less")

wilcox.test(allGLMStats$AUC_2D, allGLMStats$AUC_3D, 
            paired=TRUE, alternative = "less")

wilcox.test(allGLMStats$AIC_2D, allGLMStats$AIC_3D, 
            paired=TRUE, alternative = "less")

wilcox.test(allGLMStats$Kappa_2D, allGLMStats$Kappa_3D, 
            paired=TRUE, alternative = "greater")

wilcox.test(allGLMStats$TSS_2D, allGLMStats$TSS_3D, 
            paired=TRUE, alternative = "greater")

write.csv(allGLMStats, "~/Dropbox/MARDIGRA/Papers/Manuscript/Appendix Material/AllModelStatistics.csv", row.names = FALSE)

AUCplot <- ggplot(allGLMStats, aes(x = AUC_2D, y = AUC_3D)) +
  geom_point() + 
  scale_color_gradient() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method='lm', col = "red") +
  labs(x = "2D Model AUC", y = "3D Model 3D") +
  theme_classic() +
  theme(legend.position = "none")

AICplot <- ggplot(allGLMStats, aes(x = AIC_2D, y = AIC_3D)) +
  geom_point() + 
  scale_color_gradient() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method='lm', col = "red") +
  labs(x = "2D Model AIC", y = "3D Model AIC") +
  theme_classic() +
  theme(legend.position = "none")

TSSplot <- ggplot(allGLMStats, aes(x =TSS_2D, y = TSS_3D)) +
  geom_point() + 
  scale_color_gradient() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method='lm', col = "red") +
  labs(x = "2D Model TSS", 
       y = "3D Model TSS") +
  theme_classic() +
  theme(legend.position = "none")

KappaPlot <- ggplot(allGLMStats, aes(x =Kappa_2D, y = Kappa_3D)) +
  geom_point() + 
  scale_color_gradient() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method='lm', col = "red") +
  labs(x = "2D Model Kappa", 
       y = "3D Model Kappa") +
  theme_classic() +
  theme(legend.position = "none")

# Plot 'em
pdf("~/Dropbox/MARDIGRA/Papers/Manuscript/Appendix Material/ModelStatPlots.pdf")
ggarrange(AUCplot, AICplot,
          TSSplot, KappaPlot, 
          ncol = 2, nrow = 2)
dev.off()

# For overplotting
names(glm2D) <- ogNames
names(glm3D) <- ogNames
glm2D$modelDim <- "2D"
glm3D$modelDim <- "3D"
allGLMStatsLong <- rbind(glm2D,glm3D)

AUCpres <- ggplot(allGLMStatsLong, aes(x = nPres, y = AUC, col = modelDim)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x = "Number of Presences") +
  theme_classic() +
  theme(legend.position = "none")

AICpres <- ggplot(allGLMStatsLong, aes(x = nPres, y = AIC, col = modelDim)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x = "Number of Presences") +
  theme_classic() +
  theme(legend.position = "none")

TSSpres <- ggplot(allGLMStatsLong, aes(x = nPres, y = TSS, col = modelDim)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x = "Number of Presences") +
  theme_classic() +
  theme(legend.position = "none")

KappaPres <- ggplot(allGLMStatsLong, aes(x = nPres, y = Kappa, col = modelDim)) +
  geom_point() +
  geom_smooth(method='lm') +
  labs(x = "Number of Presences") +
  theme_classic() +
  theme(legend.position = "none")

legendPlot <- ggplot(allGLMStatsLong, aes(x = AUC, y = AUC, 
                                      col = modelDim))+
  geom_point()+
  lims(x = c(0,0), y = c(0,0))+
  theme_void()+
  labs(col="Model\nDimensions")+
  theme(legend.position = c(0.5,0.5),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size =  12),
        legend.title = element_text(size = 15, face = "bold"))+
  guides(colour = guide_legend(override.aes = list(size=8)))

# Plot 'em
pdf("~/Dropbox/MARDIGRA/Papers/Manuscript/Appendix Material/ModelStatPlotsVsPresences.pdf")
ggarrange(AUCpres, AICpres,
          TSSpres, KappaPres, legendPlot,
          ncol = 3, nrow = 2)
dev.off()