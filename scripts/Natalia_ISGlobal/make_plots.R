## Make plots for Natalia
library(minfi)
library(ggplot2)

load("results/methylation/finalQC_files/2021-03-08/gset.autosomic.Rdata")

df <- data.frame(cpg = as.vector(getBeta(gset["cg16810626", ])), status =  gset$Status)
df$status <- factor(df$status, levels = c("Control", "Case"))

png("scripts/Natalia_ISGlobal/cg16810626_chd.png")
ggplot(data = df, aes(x = status, y = cpg))+
  geom_boxplot()+
  labs(x = "Congenital heart defect", y = "cg16810626 methylation")+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()
