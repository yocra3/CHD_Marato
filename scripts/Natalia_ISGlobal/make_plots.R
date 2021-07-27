## Make plots for Natalia
library(minfi)
library(ggplot2)

load("results/methylation/finalQC_files/2021-03-08/gset.autosomic.Rdata")

df <- data.frame(cpg = as.vector(getBeta(gset["cg16810626", ])), status =  gset$Status)
df$status <- factor(df$status, levels = c("Control", "Case"))

plotA <- ggplot(data = df, aes(x = status, y = cpg))+
  geom_boxplot(fill="#BEBEBE") +
  labs(x="Congenital Heart Defect", y="cg16810626 methylation")+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        strip.text.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  annotate("rect", xmin = 0.5, xmax = 1.1, ymin = 0.5, ymax = 0.65, fill="white", colour="black") +
  annotate("text", y=0.60, x=0.8, label = "P = 4.4e-05")+
  annotate("text", y=0.54, x=0.8, label = "beta == 0.153", parse=T)


png("scripts/Natalia_ISGlobal/cg16810626_chd.v2.png", res = 1200, width = 5, height = 5, units = "in")
par(mar=c(5,5,2,2))
print(plotA)
dev.off()