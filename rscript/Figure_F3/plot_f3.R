# Load libraries
suppressWarnings(library("gridExtra"))
suppressWarnings(library("ggplot2"))
#install.packages('Hmisc')
library(Hmisc)

# Define colour palette for meta-populations
palette <- list("#a6cee3"="African Dogs",
                "#cab2d6"="American Dogs",
                "#b15928"="American Wolf", 
                "#4d4d4d"="Ancient Wolf",
                "#fb9a99"="Arctic Dogs", 
                "#fdbf6f"="Asian Dogs",
                "#1f78b4"="Coyotes", 
                "#b2df8a"="CTVT",
                "#ff7f00"="Dingo",
                "#e31a1c"="East Asian Dogs",
                "#003c30"="Eurasian Wolf",
                "#33a02c"="European Dogs",
                "#4d4d4d"="Outgroup", 
                "#6a3d9a"="Pre-Columbian Dogs",
                "#FFFF00"="Spitz Dogs")

# Define function to allow us to pull out colours for subsets of meta populations
f <- function(x) as.character(names(palette[match(x, palette)]))

# High quality DPC samples
dpc.highqual<-c("AL3223","AL3194")


  pop_names <-read.csv('pop_names.csv', header=TRUE)
  fam <- read.table('dpc_genome.fam', header=TRUE)
  pops <- merge(pop_names, fam[,c(1,2)], by.x='Code', by.y='FID')
  # need to replace using levels
  levels(pops$IID)[46] <-"TAI"
  
  ##################
  # GENOME
  ##################
  
  # Read in F3 stats 
  f3.genome <- read.table('f3_genome.txt', header=FALSE)
  colnames(f3.genome)=c("PopB", "PopC", "F3")  
  f3.genome.pops<- merge(f3.genome, pops, by.x='PopC', by.y='IID')
  
  # Fix target population
  for(target in dpc.highqual){
    assign(paste0("f3.genome.",target), f3.genome.pops[which(f3.genome.pops$PopB==target),])
    f3.genome.AL3194 <- f3.genome.AL3194[!grepl('AL3223', f3.genome.AL3194$PopC),]
    
    # Make new palette (NA wolf is excluded and this messes up the colours)
    f3.palette <- f(c(as.character(sort(unique(get(paste0("f3.genome.",target))$Type.Name)))))
    
    assign(paste0("f3.plot.",target), ggplot(get(paste0("f3.genome.",target)), aes(reorder(PopC, F3), F3, colour=Type.Name)) +
             xlab("Population") + 
             ylab("f3") +
             coord_flip() +
             geom_point(size=3.5) +
             geom_point(shape=1, size=3.5, colour="black") +
             scale_colour_manual(values=f3.palette) +
             labs(title=paste0("f3(Outgroup;",target,",X)")) +
             theme_bw() +
             theme(legend.title=element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),text = element_text(size=14)))      
    # To add standard error bars
    #geom_linerange(aes(ymin=F3-StdErr, ymax=F3+StdErr), col="black",lwd=2) +
    
    # Hacky. Make gradient palette to match F3 heatmap
    f3.gradient.palette <- c(rep("#3800C6",50),rep("#0000FF",50),rep("#7984d6",50),rep("#dbdb51", 30),rep("#FFFF00",30),rep("#FF0000",20))
    
    assign(paste0("f3.plot.gradient.",target), ggplot(get(paste0("f3.genome.",target)), aes(reorder(PopC, F3), F3, colour=F3)) +
             xlab("Population") + 
             ylab("f3") +
             coord_flip() +
             geom_point(size=3.5) +
             geom_point(shape=1, size=3.5, colour="black") +
             scale_colour_gradientn(colors=f3.gradient.palette) +
             labs(title=paste0("f3(Outgroup;",target,",X)")) +
             theme_bw() +
             theme(legend.title=element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),text = element_text(size=14)))      
    # To add standard error bars
    #geom_linerange(aes(ymin=F3-StdErr, ymax=F3+StdErr), col="black",lwd=2) +
  }

pdf('f3_genome_meta.pdf', useDingbats=FALSE, width=14, height=12)
grid.arrange(f3.plot.AL3194 + theme(legend.position="none"), f3.plot.AL3223 ,ncol=2, widths=c(2.05, 2.6))    
dev.off()

pdf('f3_genome_gradient.pdf', useDingbats=FALSE, width=10, height=12)
grid.arrange(f3.plot.gradient.AL3194 + theme(legend.position="none"), f3.plot.gradient.AL3223 ,ncol=2,widths=c(2.05, 2.6))
dev.off()

##################
# SNP
##################

f3.snps <- read.table('f3_SNParray.txt', header=FALSE)
f3.pop <- read.table('f3_pop.out', header=TRUE, sep=',')
colnames(f3.snps)=c("PopC", "F3")
f3.snps.pops <- merge(f3.snps, f3.pop, by.x='PopC', by.y='Breed')

f3.plot.snp <- ggplot(f3.snps.pops, aes(reorder(PopC, F3), F3, colour=IID)) +
  xlab("Population") + 
  ylab("f3") +
  coord_flip() +
  geom_point(size=2.2) +
  geom_point(shape=1, size=2.2, colour="black") +
  scale_colour_manual(values=f3.palette) +
  theme_bw() +
  theme(legend.title=element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),text = element_text(size=10))    
# Add if you want to include standard error bars
#geom_linerange(aes(ymin=F3-StdErr, ymax=F3+StdErr), col="black",lwd=2) +

pdf('f3_snp_meta.pdf', useDingbats = FALSE, height=22, width=10)
f3.plot.snp
dev.off()

f3.plot.gradient.snp <- ggplot(f3.snps.pops, aes(reorder(PopC, F3), F3, colour=F3)) +
  xlab("Population") + 
  ylab("f3") +
  coord_flip() +
  geom_point(size=2.2) +
  geom_point(shape=1, size=2.2, colour="black") +
  scale_colour_gradientn(colors=f3.gradient.palette) +
  theme_bw() +
  theme(legend.title=element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),text = element_text(size=10))    
# Add if you want to include standard error bars
#geom_linerange(aes(ymin=F3-StdErr, ymax=F3+StdErr), col="black",lwd=2) +

pdf('f3_snp_gradient.pdf', useDingbats = FALSE, height=22, width=10)
f3.plot.gradient.snp
dev.off()
