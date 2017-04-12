# Last updated 11-04-17
#install.packages("gridExtra")
#install.packages("ggplot2")
library(gridExtra)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
#pops <- list(AfricanDogs=c("Basenji","Namibian Village Dog"), AmericanWolf=c("American Wolf"),NorthernDogs=c("")

if (length(args)==3){
  dstat1 <- read.table(args[1], header=FALSE)
  dstat2 <- read.table(args[2], header=FALSE)
  pdf(args[3], useDingbats = FALSE)
  panel1 <- ggplot(dstat1, aes(V5, V6)) +
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot(position=position_dodge(width=0.8)) +
  coord_flip() +
  ggtitle("a") +
  xlab("") + 
  ylab("D") +
  geom_hline(yintercept=0.0, colour="red", linetype="dashed", size=1 ) +
  theme_bw() + theme(panel.border =  element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  panel2 <-ggplot(dstat2, aes(V5, V6)) +
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot(position=position_dodge(width=0.8)) +
  coord_flip() +
  ggtitle("b") +
  xlab("") + 
  ylab("D") +
  geom_hline(yintercept=0.0, colour="red", linetype="dashed", size=1 ) +
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

grid.arrange(panel1, panel2, ncol=2)
dev.off()
} else if (length(args)==2){
  dstat1 <- read.table(args[1], header=FALSE)
  pdf(args[2], useDingbats = FALSE)
  panel1 <- ggplot(dstat1, aes(V5, V6)) +
    stat_boxplot(geom ='errorbar') + 
    geom_boxplot(position=position_dodge(width=0.8)) +
    coord_flip() +
    ggtitle("a") +
    xlab("") + 
    ylab("D") +
    geom_hline(yintercept=0.0, colour="red", linetype="dashed", size=1 ) +
    theme_bw() + theme(panel.border =  element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  grid.arrange(panel1, ncol=1)
  dev.off()
} else if (length(args)<2){
  print("Required arguments: dstats input and name for output pdf.")
}
