##### NMDS plot with vectors
###Data Input: virome - Relative read counts of viruses/ bacteria with cohort label
library(ggplot2)
library(vegan)
#select state
virome_env <- virome[1:3]
#remove env factor
virome1 <- virome[-c(1:3)]

#mds and envfit functions
set.seed(123)
virome.mds <- metaMDS(virome1, distance = "bray", autotransform = F)
virome.envfit <- envfit(virome.mds, virome_env, permutations = 99999) 

#Intrinsic variables: investigate the species which may be driving the site distribution pattern
#create for vectors construction: species vectors
virome.spp.fit <- envfit(virome.mds, virome1, permutations = 99999)

#save NMDS results into dataframe
site.scrs <-as.data.frame(scores(virome.mds, display = "sites"))

#add grouping variable "state" to dataframe
site.scrs <-cbind(site.scrs, State = virome_env$State)
site.scrs <-cbind(site.scrs, Age = virome_env$Age.Range)
site.scrs <-cbind(site.scrs, Race = virome_env$Race)

#save species intrinic values into dataframe
spp.scrs <-as.data.frame(scores(virome.spp.fit, display = "vectors"))
#add species names to dataframe
spp.scrs <-cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs <-cbind(spp.scrs, pval = virome.spp.fit$vectors$pvals)
sig.spp.scrs <- subset(spp.scrs, pval <= 0.001) 

#plotting of ordination plot
site.scrs$Race <- factor(site.scrs$Race, levels = c("Chinese", "Malay", "Indian", "Others", "Not Stated"))
Accent=c("#7fc97f","#beaed4", "#fdc086","#386cb0","#f0027f","#bf5b17")
nmds.plot.virome <-ggplot(site.scrs, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = factor(State), shape = factor(State)), size = 2)+
  coord_fixed()+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "State", shape = "State") +
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))+
  scale_color_manual(values = Accent) 
nmds.plot.virome #plot

#plot significant species
nmds.plot.virome + 
  geom_segment(data = sig.spp.scrs, aes(x=0, xend=NMDS1, y = 0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = "#E96245", lwd = 0.5, alpha = 0.8) +
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)

#ANOSIM test 
library(vegan)
com = virome[,4:ncol(virome)]
m_com = as.matrix(com)
ano = anosim(m_com, virome$State, distance = "bray", permutations = 9999)
ano  

####Boxplot
##Data Input: ttmydf- Relative read counts of viruses/ bacteria with cohort labelled
library(reshape2) 
ttmydf[ttmydf == 0] <- NA 
mdatanona=ttmydf[, -which(colMeans(is.na(df)) > 0.67)] ##Remove viruses with more than 67% missing values
amdatanona <- melt(df, id=c("ID","State")) 
library(ggplot2) 
data <- amdatanona 
data$state <- as.factor(data$state) 
Accent=c("#7fc97f","#beaed4", "#fdc086")
Dark2=c("#1b9e77","#7570b3","#d95f02") 
ggplot(data, aes(x = state, y = value)) + theme_minimal()+ theme(strip.text = element_text(face = "italic"))+ 
  facet_wrap(~variable, scales = "free",ncol = 6)+  
  xlab("state") + ylab("Log 10 Relative Abundance")+ 
  geom_boxplot(aes(fill = state), alpha = 0.2 , position = position_dodge2(preserve = "total"), outlier.shape = NA) + 
  geom_jitter(aes(color = state))+ 
  scale_fill_manual(values=Accent)+ 
  scale_color_manual(values = Dark2)

#### Interbacterial correlation plot
##Data input: tmydf – Relative read counts of virus with cohort labelled 
tmydf[tmydf ==0] <- NA 
mdatanona=tmydff[, -which(colMeans(is.na(tmydff)) > 0.67)] ##remove bacteria with more than 67% missing value
amdatanona <- melt(mdatanona, id=c("ID","state")) 

##subset of dataframe based on cohort
healthy=subset(amdatanona, State=='Healthy') 
healthy <- dcast(healthy, ID ~ variable, value.var='value') 
ahfmd=subset(amdatanona, State=='Asymptomatic') 
ahfmd <- dcast(ahfmd, ID ~ variable, value.var='value') 
shfmd=subset(amdatanona, State=='Symptomatic') 
shfmd <- dcast(shfmd, ID ~ variable, value.var='value') 

##Remove columns with only zeros
healthy[healthy == 0] <- NA 
healthy=healthy[, -which(colMeans(is.na(healthy)) > 0.05)] 
ahfmd[ahfmd == 0] <- NA 
ahfmd=ahfmd[, -which(colMeans(is.na(ahfmd)) > 0.05)] 
shfmd[shfmd == 0] <- NA 
shfmd=shfmd[, -which(colMeans(is.na(shfmd)) > 0.05)] 
dhealthy <- healthy[c(-1)] 
dhealthy[is.na(dhealthy)] <- 0 
dahfmd <- ahfmd[c(-1)] 
dahfmd[is.na(dahfmd)] <- 0 
dshfmd <- shfmd[c(-1)] 
dshfmd[is.na(dshfmd)] <- 0

##Generate spearman correlation matrix
library("Hmisc")
chealthy <- rcorr(as.matrix(dhealthy), type = c("spearman"))
cshfmd <- rcorr(as.matrix(dshfmd), type = c("spearman"))
cahfmd <- rcorr(as.matrix(dahfmd), type = c("spearman"))

##Plot to visualize spearman correlation in bacteriome community
library(corrplot)
corrplot(chealthy$r, type="upper", order="hclust", tl.col = "black", tl.srt = 90,tl.cex = 0.3,
         p.mat = chealthy$P, sig.level = 0.05, insig = "blank") ##Healthy interbacterial plot 
corrplot(cshfmd$r, type="upper", order="hclust", tl.col = "black", tl.srt = 90,tl.cex = 0.3, 
         p.mat = cshfmd$P, sig.level = 0.05, insig = "blank") ##Symptomatic interbacterial plot 
corrplot(cahfmd$r, type="upper", order="hclust", tl.col = "black", tl.srt = 90,tl.cex = 0.3, 
         p.mat = cahfmd$P, sig.level = 0.05, insig = "blank") ##Asymptomatic interbacterial plot

#### Correlation chart between enterovirus and bacterial community
##Data Input: dataAS – Relative read counts of enterovirus and bacteria 

##Generate spearman correlation matrix 
library("Hmisc") 
res2 <- rcorr(as.matrix(dataAS[-c(1,2)]), type = c("spearman")) 
flattenCorrMatrix <- function(cormat, pmat) { 
  ut <- upper.tri(cormat) 
  data.frame( 
    row = rownames(cormat)[row(cormat)[ut]], 
    column = rownames(cormat)[col(cormat)[ut]], 
    cor =(cormat)[ut], 
    p = pmat[ut] 
  ) 
} 
res3=flattenCorrMatrix(res2$r, res2$P) 
res4=subset(res3, row == “Enterovirus A”) 
res45=subset(res4, p<0.05) ##Select for significant correlation between enterovirus and bacteria 

##Plot only significant correlated enterovirus and bacteria 
library(PerformanceAnalytics) 
sigdataAS = dataAS[c(1,3,19,46,47,62,65,68)] 
chart.Correlation(sigdataAS, histogram=TRUE, pch=19, method=“spearman”)