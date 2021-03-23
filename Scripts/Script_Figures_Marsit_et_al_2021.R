##R version 4.0.4  was used to run the following scripts

#####Libraries used

#flowCore_1.42.3
library("flowCore")
##flowViz_1.40.0
library("flowViz")
##ggplot2_3.3.1
library(ggplot2)
##lattice_0.20-35
library(lattice)
##gridExtra_2.3
library(gridExtra)
##plyr_1.8.4
library(plyr)
##reshape2_1.4.3
library(reshape2)
##lemon_0.4.1
library(lemon)
##tidyverse_1.3.0
library(tidyverse)
##dplyr_1.0.0
library(dplyr)
##tidyr_1.1.0
library(tidyr)
##data.table_1.12.0
library(data.table)
##magrittr_1.5
library(magrittr)
##stringr_1.4.0
library(stringr)
##UpSetR_1.4.0
library(UpSetR)
##grid_3.4.1
library(grid)
##ggpubr_0.2.2
library(ggpubr)
##gdata_2.18.0
library(gdata)
#filesstrings_3.0.0
library(filesstrings)
##cowplot_0.9.3
library(cowplot)
#ggbeeswarm_0.6.0
library("ggbeeswarm")

#Below are the scripts used for figures. All the data tables are in the Source_Data file, transform first all tables in .txt files
#Put the path to the folder where all the .txt are
setwd("")

##################################################################### Figure 1 B ploidy kinetics #######################################################

#read tables with fluorescence data for each cross and order them for the figure

data_P_all = read.table("tab_Fig_1_B.txt",header=T,sep="\t")

data_P_all $Generations=as.numeric(data_P_all $Generations)
data_P_all $cross_rep=as.factor(data_P_all $cross_rep)
data_P_all $Total_4n=as.numeric(data_P_all $Total_4n)

data_P_all_1= filter (data_P_all, Hyb == "VL_B" | Hyb == "VL_A" | Hyb == "VL_S")
data_P_all_1 $Hyb_f = factor(data_P_all_1 $Hyb, levels=c("VL_B", "VL_A", "VL_S"))
data_P_all_1 $type_f = factor(data_P_all_1 $type, levels=c("Parents 1"))

data_P_all_2= filter (data_P_all, Hyb == "VL_C" | Hyb == "VL_A" | Hyb == "VL_S")
data_P_all_2 $Hyb_f = factor(data_P_all_2 $Hyb, levels=c("VL_C", "VL_A", "VL_S"))
data_P_all_2 $type_f = factor(data_P_all_2 $type, levels=c("Parents 2"))

data_P_all_3= filter (data_P_all, Hyb == "L" | Hyb == "M" | Hyb == "H")
data_P_all_3 $Hyb_f = factor(data_P_all_3 $Hyb, levels=c("L", "M", "H"))
data_P_all_3 $type_f = factor(data_P_all_3 $type, levels=c("Hybrids"))


blue_GC= rgb(0,0,(139/255))
red_GC= rgb((238/255),0,0)
green_GC= rgb(0,(139/255),0)
black_GC= rgb(0,0,0)

pdf(("Fig_1_B.pdf"), height = 10, width =9 )
fig3B <-
  ggplot(data_P_all_1, aes(x = as.numeric(Generations), y = as.numeric(Total_4n), group=cross)) + 
  geom_line (aes(linetype=cross_rep)) +
  scale_linetype_manual(values=c("solid", "dashed" , "dotted"))+
  theme(legend.position="top")+
  ylab("")+ theme(axis.text=element_text(size=14))+
  xlab("")+ theme(axis.text=element_text(size=14))+
  facet_wrap(~ Hyb_f, ncol=3)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
  theme(axis.text.x=element_text(size=12, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black"))+
  ylim(0, 7)

g <- ggplot_gtable(ggplot_build(fig3B))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c(red_GC)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid::grid.draw(g)


fig3C <-
  ggplot(data_P_all_2, aes(x = as.numeric(Generations), y = as.numeric(Total_4n), group=cross)) + 
  geom_line (aes(linetype=cross_rep)) +
  scale_linetype_manual(values=c("solid", "dashed" , "dotted"))+
  theme(legend.position="top")+
  ylab("Number of tetraploids")+ theme(axis.text=element_text(size=14))+
  xlab("")+ theme(axis.text=element_text(size=14))+
  facet_wrap(~ Hyb_f, ncol=3)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
  theme(axis.text.x=element_text(size=12, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black"))+
  ylim(0, 7)

h <- ggplot_gtable(ggplot_build(fig3C))
stripr <- which(grepl('strip-t', h$layout$name))
fills <- c("deepskyblue3","yellow4" ,"maroon4")

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', h$grobs[[i]]$grobs[[1]]$childrenOrder))
  h$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid::grid.draw(h)


fig3D <-
  ggplot(data_P_all_3, aes(x = as.numeric(Generations), y = as.numeric(Total_4n), group=cross)) + 
  geom_line (aes(linetype=cross_rep)) +
  scale_linetype_manual(values=c("solid", "dashed" , "dotted"))+
  theme(legend.position="top")+
  ylab("")+ theme(axis.text=element_text(size=14))+
  xlab("Generations")+ theme(axis.text=element_text(size=14))+
  facet_wrap(~ Hyb_f, ncol=3)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
  theme(axis.text.x=element_text(size=12, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black"))+
  ylim(0, 7)

m <- ggplot_gtable(ggplot_build(fig3D))
stripr <- which(grepl('strip-t', m$layout$name))
fills <- c(blue_GC, green_GC, black_GC)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', m$grobs[[i]]$grobs[[1]]$childrenOrder))
  m$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid::grid.draw(m)

cowplot::plot_grid(g, h, m, nrow=3)

dev.off()


############################################################# Figure 1 C ploidy frequency #######################################################################################

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)

Ploid_freq = read.table("tab_Fig_1_C.txt",header=T,sep="\t")

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

color<- c(red, red, red, "deepskyblue3", "deepskyblue3", "deepskyblue3", "maroon4", blue ,blue, green, green, black ,black)
Ploid_freq1 =filter(Ploid_freq, rep != "VL1")

pdf("Fig_1_C.pdf", width=16, height= 7)
ggplot(Ploid_freq, aes(x = (Hyb), y = as.numeric(WGD_freq), fill=Hyb)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=0.75)+ scale_x_discrete(limit = c("VL_S1", "VL_S2", "VL_A", "VL_C1","VL_C2", "VL_C3", "VL_B1","VL_B2", "VL_B3", "L1", "L2", "M1", "M2", "H1", "H2")) +
  scale_fill_manual(values=c("VL_B1" = red, "VL_B2" =red, "VL_B3" =red, "VL_C1" ="deepskyblue3", "VL_C2" ="deepskyblue3", "VL_C3" ="deepskyblue3", "VL_A" = "yellow4", "VL_S1" ="maroon4", "VL_S2" ="maroon4", "L1" =blue ,"L2" =blue, "M1" =green, "M2" =green, "H1" =black , "H2" =black))+
  ylab("Frequency of tetraploids") +
  xlab("Crosses")+
  theme(axis.text.x=element_text(size=18), axis.title.x=element_text(size=24)) +
  theme(axis.text.y=element_text(size=24), axis.title.y=element_text(size=24))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()



################################################## Figure S1 #########################################################################

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)
library(lemon)

tab_P_all = read.table("tab_Fig_S1.txt",header=T,sep="\t")

#Separate the crosses according to the different panels
data_P_all_1 = filter(tab_P_all, Hyb == "VL_B1" | Hyb == "VL_C1" | Hyb == "VL_A"|  Hyb == "L1"| Hyb == "M1"| Hyb == "H1")
data_P_all_2 = filter(tab_P_all, Hyb == "VL_B2" | Hyb == "VL_C2" | Hyb == "VL_S1"|  Hyb == "L2"| Hyb == "M2"| Hyb == "H2")
data_P_all_3 = filter(tab_P_all, Hyb == "VL_B3" | Hyb == "VL_C3" | Hyb == "VL_S2"|  Hyb == "L1"| Hyb == "M1"| Hyb == "H1")

data_P_all_1 $Hyb_f = factor(data_P_all_1 $Hyb, levels=c("VL_B1", "VL_C1", "VL_A", "L1", "M1", "H1"))
data_P_all_1 $Time_f = factor(data_P_all_1 $Time, levels=c("Tini", "Tmid", "Tend"))
data_P_all_1 $Generations=as.factor(data_P_all_1 $Generations)
data_P_all_1 $normalisation=as.numeric(data_P_all_1 $normalisation)
data_P_all_2 $Hyb_f = factor(data_P_all_2 $Hyb, levels=c("VL_B2", "VL_C2", "VL_S1", "L2", "M2", "H2"))
data_P_all_2 $normalisation=as.numeric(data_P_all_2 $normalisation)
data_P_all_2 $Time_f = factor(data_P_all_2 $Time, levels=c("Tini", "Tmid", "Tend"))
data_P_all_3 $Hyb_f = factor(data_P_all_3 $Hyb, levels=c("VL_B3", "VL_C3", "VL_S2", "L1", "M1", "H1"))
data_P_all_3 $normalisation=as.numeric(data_P_all_3 $normalisation)
data_P_all_3 $Time_f = factor(data_P_all_3 $Time, levels=c("Tini", "Tmid", "Tend"))


blue_GC= rgb(0,0,(139/255))
red_GC= rgb((238/255),0,0)
green_GC= rgb(0,(139/255),0)
black_GC= rgb(0,0,0)


pdf(("Fig_S1.pdf"), height = 10, width =20)

fig3B <-
  ggplot(data_P_all_1, aes(x = Time_f, y = as.numeric(normalisation), group=ech)) + 
  geom_line(alpha=2/10) + 
  geom_point(size=2.5, fill = "black", alpha=2/10) + 
  ylab("Ploidy")+ theme(axis.text=element_text(size=18))+
  xlab("")+
  facet_rep_grid(.~Hyb_f)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
  theme(axis.text.x=element_text(size=16, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black"))
ylim(1.5,4.5)
g <- ggplot_gtable(ggplot_build(fig3B))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c(red_GC, "deepskyblue3", "yellow4", blue_GC, green_GC, black_GC)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

fig3C <-
  ggplot(data_P_all_2, aes(x = Time_f, y = as.numeric(normalisation), group=ech)) + 
  geom_line(alpha=2/10) + 
  geom_point(size=2.5, fill = "black", alpha=2/10) + 
  ylab("Ploidy")+ theme(axis.text=element_text(size=18))+
  xlab("")+
  facet_rep_grid(.~Hyb_f)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
  theme(axis.text.x=element_text(size=16, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black"))
ylim(1.5,4.5)

g1 <- ggplot_gtable(ggplot_build(fig3C))
stripr <- which(grepl('strip-t', g1$layout$name))
fills <- c(red_GC, "deepskyblue3", "maroon4", blue_GC, green_GC, black_GC)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g1)

fig3d <-
  ggplot(data_P_all_3, aes(x = Time_f, y = as.numeric(normalisation), group=ech)) + 
  geom_line(alpha=2/10) + 
  geom_point(size=2.5, fill = "black", alpha=2/10) + 
  ylab("Ploidy")+ theme(axis.text=element_text(size=18))+
  xlab("")+
  facet_rep_grid(.~Hyb_f)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
  theme(axis.text.x=element_text(size=16, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black"))
ylim(1.5,4.5)

g2 <- ggplot_gtable(ggplot_build(fig3d))
stripr <- which(grepl('strip-t', g2$layout$name))
fills <- c(red_GC, "deepskyblue3", "maroon4", "white", "white", "white")

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g2)
cowplot::plot_grid(g, g1, g2, nrow=3)
dev.off()

######################################################## Figure S2 ######################################################################

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)
library(lemon)
library(grid)
rm(list=ls())

data_P_all_1 = read.table("tab_Fig_S2.txt",header=T,sep="\t")

data_P_all_1 $cross_f = factor(data_P_all_1 $cross, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))
data_P_all_1 $Generations_f = factor(data_P_all_1 $Generations, levels=c(22, 88, 154, 220, 286, 352, 440, 460, 484, 550, 572, 616, 638, 682, 704 ,770))
data_P_all_1 $normalisation=as.numeric(data_P_all_1 $normalisation)
data_P_all_1= filter(data_P_all_1, cross_f != "NA")
data_P_all_1 $Generations=as.numeric(data_P_all_1 $Generations)

blue_GC= rgb(0,0,(139/255))
red_GC= rgb((238/255),0,0)
green_GC= rgb(0,(139/255),0)
black_GC= rgb(0,0,0)


pdf(("Fig_S2.pdf"), height = 4, width =38 )
fig3B <- 
  ggplot(data_P_all_1, aes(x = Generations, y = as.numeric(normalisation), group=strain)) + 
  geom_line(alpha=2/10) + 
  geom_point(size=2.5, fill = "black", alpha=4/10) + 
  ylab("Ploidy")+ theme(axis.text=element_text(size=18))+
  xlab("")+
  facet_rep_grid(.~cross_f)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
  theme(axis.text.x=element_text(size=12, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black"))
ylim(1.5,4.5)

g <- ggplot_gtable(ggplot_build(fig3B))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c("deepskyblue3", "deepskyblue3", "deepskyblue3", blue_GC, blue_GC, green_GC, black_GC)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g) 

dev.off()

######################################################### Figure S4 ############################################################################################

library("ggbeeswarm")
rm(list=ls())

Tab_P_heatmap1 = read.table("tab_Fig_S4.txt",header=T,sep="\t")
Tab_P_heatmap1$Lineage = factor(Tab_P_heatmap1$lineage, levels=c("Ctrl", "SpA", "SpB", "SpC"))


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)


pdf ('Fig_S4.pdf', width = 9 , height = 5)
ggplot(Tab_P_heatmap1, aes(x = Lineage, y = as.numeric(normalisation), fill = Lineage)) + geom_beeswarm(aes(color = Lineage, binwidth = 0.03), alpha=4/10) +
  scale_color_manual(values=c("grey10" , "yellow4", red, "deepskyblue3"))+
  ylab("Ploidy")+ theme(axis.title =element_text(size=18))+
  xlab("Lineages")+ theme(axis.title=element_text(size=18))+
  theme(axis.text.x=element_text(size=18, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black"))+
  theme(axis.line = element_line(colour = "black"))+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()

############################################################################################
#                                  Fertility data                                          #
############################################################################################ 

################################################### Figure 2 B ##########################################################################################################

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)
library(lemon)
library(grid)
library(ggpubr)
rm(list=ls())

#read tables with fertility data for each cross
data = read.table("tab_Fig_2B.txt",header=T,sep="\t")

data $ploidy_f = factor(data $ploidy, levels=c("2n", "4n"))
data $cross_f = factor(data $cross, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))
data $time_f = factor(data $time, levels=c("Tini", "Tend"))

data_4= filter(data, fertility != "NA")
data_4= filter(data_4, type == "tetraploid")

blue_GC= rgb(0,0,(139/255))
red_GC= rgb((238/255),0,0)
green_GC= rgb(0,(139/255),0)
black_GC= rgb(0,0,0)

pdf(("Fig_2_B.pdf"), height = 5.5, width =13 )
fig3B <-
  ggplot(data= data_4, aes(y = fertility, x = ploidy_f))+#, group=strain)) + 
  geom_point(fill = "black", alpha=4/10, size=4) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.6)+
  stat_compare_means(aes(group = ploidy_f), method=("t.test"), Paired= TRUE, label = "p.format", size =4)+
  ylab("%Fertility")+ theme(axis.text=element_text(size=18))+
  xlab("Ploidy")+ theme(axis.title.y=element_text(size=22))+
  theme(axis.title.x=element_text(size=22))+
  facet_grid(.~cross_f)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ 
  theme(axis.text.x=element_text(size=18, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black")) +
  ylim(0, 110)

g <- ggplot_gtable(ggplot_build(fig3B))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c("deepskyblue3", "deepskyblue3", "deepskyblue3", blue_GC, blue_GC, green_GC, black_GC)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g) 

dev.off()

######################################################## Figure S7 #####################################################################################################
rm(list=ls())
data = read.table("tab_Fig_S7.txt",header=T,sep="\t")

data $time_f = factor(data $time, levels=c("Tini", "Tend"))
data $cross_f = factor(data $cross, levels=c("VL_C1", "VL_C2", "VL_C3"))
data_2= filter(data, fertility != "NA")

blue_GC= rgb(0,0,(139/255))
red_GC= rgb((238/255),0,0)
green_GC= rgb(0,(139/255),0)
black_GC= rgb(0,0,0)

pdf(("Fig_S7.pdf"), height = 3, width =8 )
fig3B <-
  ggplot(data= data_2, aes(y = fertility, x = time_f))+ #, group=strain)) + 
  geom_point(fill = "black", alpha=4/10, size=4) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.6)+
  stat_compare_means(aes(group = time_f, label=..p.adj..), method = "wilcox.test", label.y = 103, Paired= TRUE, size =5) +
  ylab("%Fertility")+ theme(axis.text=element_text(size=18))+
  xlab("")+
  facet_grid(.~cross_f)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ 
  theme(axis.text.x=element_text(size=18, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black")) 

g <- ggplot_gtable(ggplot_build(fig3B))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c("deepskyblue3", "deepskyblue3", "deepskyblue3")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g) 

dev.off()

############################################################################################
#                        Growth rate measurement                                           #
############################################################################################ 

############################################## Figure S8 #################################################################################################################
rm(list=ls())

data = read.table("Fig_S8.txt",header=T,sep="\t")

data $time_f = factor(data $time, levels=c("Tini", "T2n", "Tmid", "T4n", "Tend"))
data $hyb_f = factor(data $hyb, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))
data=filter(data, time_f!="NA")
data=filter(data, hyb_f!="NA")
data $group_f = factor(data $group, levels=c("tetraploid", "diploid"))


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

pdf(("Fig_S8.pdf"), height = 4, width =22 )
#
fig3B <-
  ggplot(data, aes(x = time_f, y = as.numeric(median_max_slope), ymin=median_max_slope-mad_slope, ymax=median_max_slope+mad_slope, group=cross)) +
  geom_line(aes(linetype=group_f)) +
  geom_pointrange(aes(color = group_f)) +
  #geom_point(aes(color=fert), alpha = 0.85, size = 2.5) +
  scale_color_manual(values=c("violetred", "lightslategrey")) +
  #scale_x_discrete(labels = c("ctrl" = "ctrl", "0"="Tini", "~500"= "Tmid", "~1000"="Tend", size=18 )) +
  ylab("Growth rate(OD/h)")+ theme(axis.text=element_text(size=18))+
  xlab("")+
  facet_grid(.~hyb_f)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
  theme(axis.text.x=element_text(size=16, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black")) #+

g <- ggplot_gtable(ggplot_build(fig3B))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c("deepskyblue3", "deepskyblue3", "deepskyblue3", blue, blue, green, black)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

############################################## Figure S9_A #################################################################################################################
rm(list=ls())

data_2n_4n = read.table("tab_Fig_S9A.txt",header=T,sep="\t")

data_2n_4n=filter(data_2n_4n, time1=="T2n" |time1=="T4n")
data_2n_4n $time_f = factor(data_2n_4n $time1, levels=c("T2n", "T4n"))
data_2n_4n $hyb_f = factor(data_2n_4n $hyb, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

pdf(("Fig_S9_A.pdf"), height = 5, width =13 )
fig3B <-
  ggplot(data_2n_4n, aes(x = time_f, y = as.numeric(median_max_slope))) +#, group=hyb_f)) + 
  geom_line(aes(group = (cross)), alpha = 0.25, fill =("grey")) + 
  geom_point(fill = "black", alpha=4/10, size=4) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.62)+
  stat_compare_means(aes(group = time_f), method=("t.test"), Paired= TRUE, label = "p.format", size =4)+
  ylab("Growth rate (OD/h)")+ theme(axis.text=element_text(size=18))+
  xlab("Ploidy")+ theme(axis.title.y=element_text(size=22))+
  theme(axis.title.x=element_text(size=22))+
  facet_grid(.~hyb_f)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ 
  theme(axis.text.x=element_text(size=18, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black")) 

g <- ggplot_gtable(ggplot_build(fig3B))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c("deepskyblue3", "deepskyblue3", "deepskyblue3", blue, blue, green, black)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g) 

dev.off()


############################################## Figure S9_B #################################################################################################################

rm(list=ls())

diff_4n_tend = read.table("tab_Fig_S9B.txt",header=T,sep="\t")

diff_4n_tend $hyb_f = factor(diff_4n_tend $hyb, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)


pdf("Fig_S9B.pdf", width=14, height= 3.2)

fig3B <- ggplot(diff_4n_tend[!is.na(diff_4n_tend$group1),], aes(x = nbr_Gener, y = as.numeric(diff), group=hyb)) +
  geom_point( shape=21, colour ="black", fill= "black", size=4, alpha= 0.6) +
  geom_hline(aes(yintercept = 0), colour = 'red', linetype = "dotted", size=1)+
  facet_grid(.~hyb_f, space= "free")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+
  ylab("GR(Tend-TpostWGD) (OD/h)")+
  xlab("Generations (Tend-TpostWGD)")+
  theme(axis.title.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.text.x=element_text(size=16, colour = "black"))+ theme(axis.text.y=element_text(size=14, colour = "black"))

g <- ggplot_gtable(ggplot_build(fig3B))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c("deepskyblue3", "deepskyblue3", "deepskyblue3", blue, blue, green, black)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g) 
dev.off()

############################################## Figure 2C #################################################################################################################

rm(list=ls())

diff_2n_4n = read.table("tab_Fig_2C.txt",header=T,sep="\t")

diff_2n_4n $hyb_f = factor(diff_2n_4n $hyb, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))

t_tests2 = diff_2n_4n %>%
  group_by(hyb_f) %>% dplyr::summarise(P = t.test(diff, mu = 0)$p.value, Sig = ifelse(P < 0.05, "*", "ns"), max_diff = max(diff))

t_tests2 <- t_tests2 %>%  mutate_if(is.numeric, round, digits = 2)

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

pdf("Fig_2_C.pdf", width=7.5, height= 5)
ggplot(diff_2n_4n, aes(x = hyb_f, y = as.numeric(diff), fill= hyb_f)) +
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", blue, blue, green, black))+
  geom_point( shape=21, colour = "white", size=4, alpha= 0.75) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.62)+
  
  geom_hline(aes(yintercept = 0), colour = 'red', linetype = "dotted", size=1)+
  theme(axis.text.x=element_text(size=18, colour = "black"))+ theme(axis.text.y=element_text(size=16, colour = "black"))+
  ylab("Gain of fitness by WGD (OD/h)")+
  xlab("Cross")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+ theme_bw() + theme(panel.grid.major = element_blank(),
                                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(aes(label = P, y = 0.014), size = 4, data = t_tests2)

dev.off()


############################################## Figure 2D #################################################################################################################
rm(list=ls())

diff_4n_tend_4n = read.table("tab_Fig_2D.txt",header=T,sep="\t")

diff_4n_tend_4n $hyb_f = factor(diff_4n_tend_4n $hyb, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))

diff_4n_tend_4n1= filter (diff_4n_tend_4n, hyb != "L2" & hyb != "VL_C2") 

# Test whether each group differs from 0
t_tests = diff_4n_tend_4n1 %>% group_by(hyb_f) %>% dplyr::summarise(P = t.test(diff_Gener, mu = 0)$p.value, Sig = ifelse(P < 0.05, "*", "ns"), max_diff = max(diff_Gener))

t_tests <- t_tests %>%  mutate_if(is.numeric, round, digits = 3)

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

pdf("Fig_2_D.pdf", width=7.5, height= 5)
ggplot(diff_4n_tend_4n, aes(x = hyb_f, y = as.numeric(diff_Gener), fill= hyb_f)) +
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", blue, blue, green, black))+
  geom_point( shape=21, colour = "white", size=4, alpha= 0.75) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.62)+
  geom_hline(aes(yintercept = 0), colour = 'red', linetype = "dotted", size=1)+
  theme(axis.text.x=element_text(size=18, colour = "black"))+ theme(axis.text.y=element_text(size=16, colour = "black"))+
  ylab("Normalized gain of fitness post WGD (OD/h)")+
  xlab("Cross")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(aes(label = P, y = 1.12e-04), size = 4, data = t_tests)

dev.off()



############################################################################################
#                                   Genomic data                                           #
############################################################################################ 

####################################################### Figure 3 A  ##################################################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)
library(lattice)
library(gridExtra)
library(grid)
library(cowplot)
library(lemon)
library(ggpubr)

rm(list=ls())

Aneup = read.table("tab_Fig_3A.txt",header=T,sep="\t")

Aneup$cross_f <- factor(Aneup$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
Aneup$Hybb_f <- factor(Aneup$Hybb, levels = c("VL_C", "VL_B", "L","M", "H"))

Aneup= filter(Aneup, ploidy == "2n")

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)


pdf(file = paste0("Fig_3_A.pdf"), height = 6, width = 9)

my_comparisons <- list(c("VL_B", "VL_C"), c("VL_B", "L"), c("VL_B", "M"), c("VL_B", "H"),c("VL_C", "L"), c("VL_C", "M"), c("VL_C", "H"), c("L", "M"), c("L", "H"), c("M", "H"))
ggplot(Aneup[!is.na(Aneup$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y = 12.5, size= 4)+      # Add global p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

dev.off()



####################################################### Figure S10 A ##################################################################

rm(list=ls())

Aneup_gain = read.table("tab_Fig_S10A.txt",header=T,sep="\t")

Aneup_gain$cross_f <- factor(Aneup_gain$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
Aneup_gain$Hybb_f <- factor(Aneup_gain$Hybb, levels = c("VL_C", "VL_B", "L","M", "H"))

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf(file = paste0("Fig_S10_A.pdf"), height = 6, width = 9)

my_comparisons <- list(c("VL_B", "VL_C"), c("VL_C", "M"), c("VL_C", "H"))
ggplot(Aneup_gain[!is.na(Aneup_gain$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y = 7.5, size= 4)+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

dev.off()


####################################################### Figure S10 B ##################################################################

rm(list=ls())

Aneup_loss = read.table("tab_Fig_S10B.txt",header=T,sep="\t")

Aneup_loss$cross_f <- factor(Aneup_loss$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
Aneup_loss$Hybb_f <- factor(Aneup_loss$Hybb, levels = c("VL_C", "VL_B", "L","M", "H"))

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf(file = paste0("Fig_S10_B.pdf"), height = 6, width = 9)

my_comparisons <- list(c("VL_B", "VL_C"), c("VL_B", "L"), c("VL_B", "M"), c("VL_B", "H"), c("VL_C", "M"), c("VL_C", "H"), c("L", "H"))
ggplot(Aneup_loss[!is.na(Aneup_loss$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y =11, size= 4)+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

dev.off()



####################################################### Figure S11 A  ##################################################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)
library(lattice)
library(gridExtra)
library(grid)
library(cowplot)
library(lemon)
library(ggpubr)

rm(list=ls())

Aneup = read.table("tab_Fig_S11A.txt",header=T,sep="\t")

Aneup$cross_f <- factor(Aneup$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
Aneup$Hybb_f <- factor(Aneup$Hybb, levels = c("VL_C", "VL_B", "L","M", "H"))

Aneup= filter(Aneup, ploidy == "2n")

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)


pdf(file = paste0("Fig_S11_A.pdf"), height = 6, width = 9)

my_comparisons <- list(c("VL_B", "VL_C"), c("VL_B", "L"), c("VL_B", "M"), c("VL_B", "H"),c("VL_C", "L"), c("VL_C", "M"), c("VL_C", "H"), c("L", "M"), c("L", "H"), c("M", "H"))
ggplot(Aneup[!is.na(Aneup$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y = 12.5, size= 4)+      # Add global p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

dev.off()


####################################################### Figure S11 B ##################################################################

rm(list=ls())

Aneup_gain = read.table("tab_Fig_S11B.txt",header=T,sep="\t")

Aneup_gain$cross_f <- factor(Aneup_gain$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
Aneup_gain$Hybb_f <- factor(Aneup_gain$Hybb, levels = c("VL_C", "VL_B", "L","M", "H"))

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf(file = paste0("Fig_S11_B.pdf"), height = 6, width = 9)

my_comparisons <- list(c("VL_B", "VL_C"), c("VL_C", "M"), c("VL_C", "H"))
ggplot(Aneup_gain[!is.na(Aneup_gain$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y = 7.5, size= 4)+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

dev.off()


####################################################### Figure S11 C ##################################################################

rm(list=ls())

Aneup_loss = read.table("tab_Fig_S11C.txt",header=T,sep="\t")

Aneup_loss$cross_f <- factor(Aneup_loss$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
Aneup_loss$Hybb_f <- factor(Aneup_loss$Hybb, levels = c("VL_C", "VL_B", "L","M", "H"))

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf(file = paste0("Fig_S11_C.pdf"), height = 6, width = 9)

my_comparisons <- list(c("VL_B", "VL_C"), c("VL_B", "L"), c("VL_B", "M"), c("VL_B", "H"), c("VL_C", "M"), c("VL_C", "H"), c("L", "H"))
ggplot(Aneup_loss[!is.na(Aneup_loss$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y =11, size= 4)+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)+
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

dev.off()

############################################################################################### Figure 4 A ####################################

rm(list=ls())

Aneup = read.table("tab_Fig_4A.txt",header=T,sep="\t")

Aneup$cross_f <- factor(Aneup$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "L1", "L2","M1", "H2"))
Aneup$Hybb_f <- factor(Aneup$Hybb, levels = c("VL_C", "L","M", "H"))

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf(file = paste0("Fig_4_A.pdf"), height = 6, width = 9)

my_comparisons <- list( c("2n", "4n"))

ggplot(Aneup[!is.na(Aneup$cross_f),], aes(x = ploidy, y = as.numeric(Aneup_tot_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3",  blue, blue, green, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4) +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

dev.off()

########################################################### Figure S13 A ########################################################################

rm(list=ls())

Aneup = read.table("tab_Fig_S13A.txt",header=T,sep="\t")

Aneup$cross_f <- factor(Aneup$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "L1", "L2","M1", "H2"))
Aneup$Hybb_f <- factor(Aneup$Hybb, levels = c("VL_C", "L","M", "H"))

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf(file = paste0("Fig_S13_A.pdf"), height = 6, width = 9)

my_comparisons <- list( c("2n", "4n"))

ggplot(Aneup[!is.na(Aneup$cross_f),], aes(x = ploidy, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3",  blue, blue, green, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4) +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()

########################################################### Figure S13 B ########################################################################

rm(list=ls())

Aneup = read.table("tab_Fig_S13B.txt",header=T,sep="\t")

Aneup$cross_f <- factor(Aneup$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "L1", "L2","M1", "H2"))
Aneup$Hybb_f <- factor(Aneup$Hybb, levels = c("VL_C", "L","M", "H"))

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf(file = paste0("Fig_S13_B.pdf"), height = 6, width = 9)

my_comparisons <- list( c("2n", "4n"))

ggplot(Aneup[!is.na(Aneup$cross_f),], aes(x = ploidy, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3",  blue, blue, green, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4) +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()


######################################################################################## Figure 3 B ###################################################################################

rm(list=ls())

tab_freq_aneup_2n_gain = read.table("tab_Fig_3Bgain.txt",header=T,sep="\t")
tab_freq_aneup_2n_loss = read.table("tab_Fig_3Bloss.txt",header=T,sep="\t")

tab_freq_aneup_2n_gain$Hybb_f = factor(tab_freq_aneup_2n_gain$Hyb_f, levels=c("H2_2n", "H1_2n", "M2_2n", "M1_2n", "L2_2n", "L1_2n", "VL_B2_2n", "VL_B1_2n", "VL_C3_2n", "VL_C2_2n", "VL_C1_2n"))
tab_freq_aneup_2n_loss$Hybb_f = factor(tab_freq_aneup_2n_loss$Hyb_f, levels=c("H2_2n", "H1_2n", "M2_2n", "M1_2n", "L2_2n", "L1_2n", "VL_B2_2n", "VL_B1_2n", "VL_C3_2n", "VL_C2_2n", "VL_C1_2n"))


pdf(file = paste0("Fig_3_B.pdf"), height = 8, width = 20)

g_2n_g <- ggplot(tab_freq_aneup_2n_gain, aes(x = Chr , y = Hybb_f)) +
  geom_raster(aes (fill=freq_aneup_chr)) +  scale_x_discrete(expand = c(0, 0))+ scale_y_discrete(expand = c(0, 0))+
  scale_fill_gradient2(low = "white", high = "midnightblue", limits=c(0, 0.7)) + 
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25)) +
  ylab("") +xlab("")+ theme(legend.position = "bottom")
coord_equal() 

g_2n_l <- ggplot(tab_freq_aneup_2n_loss, aes(x = Chr , y = Hybb_f)) +
  geom_raster(aes (fill=freq_aneup_chr)) +  scale_x_discrete(expand = c(0, 0))+ scale_y_discrete(expand = c(0, 0))+
  scale_fill_gradient2(low = "white", high = "darkred", limits=c(0, 0.7)) + 
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25)) +
  ylab("") +xlab("")+ theme(legend.position = "bottom")
coord_equal() 
cowplot::plot_grid(g_2n_g, g_2n_l, ncol=2, rel_heights = c(300, 300), rel_widths = c(600, 600)) 

dev.off()

######################################################################################### Figure 4 B ##################################################################################

rm(list=ls())

tab_freq_aneup_4n_gain = read.table("tab_Fig_4Bgain.txt",header=T,sep="\t")
tab_freq_aneup_4n_loss = read.table("tab_Fig_4Bloss.txt",header=T,sep="\t")

tab_freq_aneup_4n_gain$Hybb_f = factor(tab_freq_aneup_4n_gain$Hyb_f, levels=c( "H2_4n", "M1_4n", "L2_4n", "L1_4n", "VL_C3_4n", "VL_C2_4n", "VL_C1_4n"))
tab_freq_aneup_4n_loss$Hybb_f = factor(tab_freq_aneup_4n_loss$Hyb_f, levels=c( "H2_4n", "M1_4n", "L2_4n", "L1_4n", "VL_C3_4n", "VL_C2_4n", "VL_C1_4n"))

pdf(file = paste0("Fig_4_B.pdf"), height = 8, width = 20)
g_4n_g <- ggplot(tab_freq_aneup_4n_gain, aes(x = Chr , y = Hybb_f)) +
  geom_raster(aes (fill=freq_aneup_chr)) +  scale_x_discrete(expand = c(0, 0))+ scale_y_discrete(expand = c(0, 0))+
  scale_fill_gradient2(low = "white", high = "midnightblue", mid = "white", limits=c(0, 1)) + 
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25)) +
  ylab("") +xlab("")+
  coord_equal() 

g_4n_l <- ggplot(tab_freq_aneup_4n_loss, aes(x = Chr , y = Hybb_f)) +
  geom_raster(aes (fill=freq_aneup_chr)) +  scale_x_discrete(expand = c(0, 0))+ scale_y_discrete(expand = c(0, 0))+
  scale_fill_gradient2(low = "white", high = "darkred", mid = "white", limits=c(0, 1)) + 
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25)) +
  ylab("") +xlab("")+
  coord_equal() 
cowplot::plot_grid(g_4n_g, g_4n_l, ncol=2, rel_heights = c(300, 300), rel_widths = c(600, 600)) 

dev.off()


############################################################## Figure 3 C #############################################################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)
library(lattice)
library(gridExtra)
library(grid)
library(cowplot)
#library(plyr)
library(lemon)
library(ggpubr)
rm(list=ls())


LOH1 = read.table("tab_Fig_3C_D.txt",header=T,sep="\t")

LOH1$cross_f <- factor(LOH1$cross, levels = c("VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
LOH1$Hyb_f <- factor(LOH1$Hyb, levels = c("VL_B", "L","M", "H"))

LOH1= filter(LOH1, ploidy == "2n")

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

pdf(file = paste0("Fig_3_C.pdf"), height = 10, width = 10)

my_comparisons <- list(c("VL_B", "L"), c("VL_B", "M"), c("VL_B", "H"), c("L", "M"), c("L", "H"))
ggplot(LOH1[!is.na(LOH1$cross_f),], aes(x = Hyb_f, y = as.numeric(n), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.2)+
  scale_fill_manual(values=c(red, red,  blue, blue, green, green, black, black))+
  facet_grid(position~., space= "free")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'black', size = 20))+
  ylab("LOH / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hyb_f),label.y = 20, size= 4)+      # Add global p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)

dev.off()


############################################################### Figure 3 D ##################################################################################
rm(list=ls())

LOH1 = read.table("tab_Fig_3C_D.txt",header=T,sep="\t")

LOH1$cross_f <- factor(LOH1$cross, levels = c("VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
LOH1$Hyb_f <- factor(LOH1$Hyb, levels = c("VL_B", "L","M", "H"))

LOH1= filter(LOH1, ploidy == "2n")


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

pdf(file = paste0("Fig_3_D.pdf"), height = 10, width = 10)

my_comparisons <- list(c("VL_B", "L"), c("VL_B", "M"), c("VL_B", "H"), c("L", "M"), c("L", "H"), c("M", "H"))

ggplot(LOH1[!is.na(LOH1$cross_f),], aes(x = Hyb_f, y = as.numeric(logsize), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.35)+
  scale_fill_manual(values=c(red, red,  blue, blue, green, green, black, black))+
  facet_grid(position~., space= "free")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'black', size = 20))+
  ylab("log2(LOH size)")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hyb_f),label.y = 18, size= 6)+      # Add global p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)+     
  stat_compare_means( method = "wilcox.test",  label = "p.format", size= 4, label.y=10.3)   
dev.off()


####################################################################################################### Figure 4 C ############################################

rm(list=ls())

LOH1 = read.table("tab_Fig_4C_D.txt",header=T,sep="\t")


LOH1$cross_f <- factor(LOH1$cross, levels = c("L1", "L2","M1", "H2"))
LOH1$Hyb_f <- factor(LOH1$Hyb, levels = c("L","M", "H"))

LOH.summary <- LOH1 %>% group_by(cross_f, Hyb, ploidy) %>% dplyr::summarise(mean=mean(n), sd = sd(n, na.rm = TRUE)) 
LOH1= filter(LOH1, cross != "VL_B1"& cross != "VL_B2"&  cross != "M2"&  cross != "H1")


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

pdf(file = paste0("Fig_4_C.pdf"), height = 7, width = 7)

my_comparisons <- list( c("2n", "4n"))
ggplot(LOH1[!is.na(LOH1$cross_f),], aes(x = ploidy, y = as.numeric(n), fill=cross_f)) +
  geom_boxplot(alpha = 0.4, color="black", size=0.25)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.35)+
  scale_fill_manual(values=c(blue, blue, green, black))+
  facet_grid(position~., space= "free")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'black', size = 20))+
  ylab("LOH / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)

dev.off()

####################################################################################################### Figure 4 D ############################################

pdf(file = paste0("Fig_4_D.pdf"), height = 7, width = 7)

my_comparisons <- list( c("2n", "4n"))
ggplot(LOH1[!is.na(LOH1$cross_f),], aes(x = ploidy, y = as.numeric(logsize), fill=cross_f)) +
  geom_boxplot(alpha = 0.4, color="black", size=0.25)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.35)+
  scale_fill_manual(values=c(blue, blue, green, black))+
  facet_grid(position~., space= "free")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'black', size = 20))+
  ylab("log2(LOH size)")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)     

dev.off()


####################################################################################################### Figure S14 ############################################

rm(list=ls())

LOH1 = read.table("tab_Fig_S14.txt",header=T,sep="\t")


LOH1$cross_f <- factor(LOH1$cross, levels = c("L1", "L2"))
LOH1$Hyb_f <- factor(LOH1$Hyb, levels = c("L","M", "H"))

LOH.summary <- LOH1 %>% group_by(cross_f, Hyb, ploidy) %>% dplyr::summarise(mean=mean(n), sd = sd(n, na.rm = TRUE)) 

LOH1= filter(LOH1, Hyb == "L")
LOH1= filter(LOH1, ploidy == "2n"| ploidy == "3n")


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)


pdf(file = paste0("Fig_S14.pdf"), height = 7, width = 7)

ggplot(LOH1[!is.na(LOH1$ploidy),], aes(x = cross_f, y = as.numeric(n), fill=ploidy)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c(blue, blue, blue, blue))+
  facet_grid(position~., space= "free")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'black', size = 20))+
  ylab("LOH / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means( method = "wilcox.test",  label = "p.format", size= 4, label.y=10.5)  

ggplot(LOH1[!is.na(LOH1$ploidy),], aes(x = cross_f, y = as.numeric(logsize), fill=ploidy)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c(blue, blue, blue, blue))+
  facet_grid(position~., space= "free")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'black', size = 20))+
  ylab("log2(LOH size)")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means( method = "wilcox.test",  label = "p.format", size= 4, label.y=10.5)  

dev.off()






