
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


############################################################################################
#                        Flow cytometry data analysis                                      #
############################################################################################ 

##########################################
# to transform a .fcs file to a .csv file#
##########################################

#source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")
biocLite("flowViz")
library("flowCore")
library("flowViz")
library(ggplot2)
library(lattice)
library(gridExtra)
library(dplyr)

#function to rename strains with their real number in the plate

num2coords = function (file_number) {
  ligne = (file_number-1) %/% 12 + 1
  colonne= (file_number-1) %% 12 +1 
  return(list (x= colonne, y = ligne))
}

coords2strain = function(list) {
  strain = ((list$x - 1)*8) + list$y
  return (strain)
}

####Transform FCS files in to CSV files

setwd("")

##### Control 1 diploid
d = read.FCS("controls/2017-07-14_ploidy controls-v2-1.fcs")

df = as.data.frame(exprs(d)[,1:ncol(d)])
ech = rep("Diplo", nrow(df))
ctrl1 = cbind(ech, df)
write.csv(ctrl1,"controls/Control1.csv", row.names = FALSE, quote = F)
csv_ctrl1 = read.table("controls/Control1.csv", header = T, sep = ",")


##### Control 2 haploid
d = read.FCS("controls/2017-07-14_ploidy controls-v2-2.fcs")

df = as.data.frame(exprs(d)[,1:ncol(d)])
ech = rep("Haplo", nrow(df))
ctrl2 = cbind(ech, df)
write.csv(ctrl2,"controls/Control2.csv", row.names = FALSE, quote = F)
csv_ctrl2 = read.table("controls/Control2.csv", header = T, sep = ",")

# Tranform samples files from FCS to CSV files with the good strain number
########################################################################

setwd("")

file = dir(pattern="fcs")
file

for(i in 1:length(file)){
  
  d = read.FCS(file[i])
  df = as.data.frame(exprs(d)[,1:ncol(d)])
  
  file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% tail(n=1) %>% as.numeric()
  strain_number= coords2strain(num2coords(file_number))
  
  ech1 = cbind(strain_number, df)
  csv_name = paste0(strsplit(file[i], ".", fixed = T)[[1]][1], "_strain_",strain_number, ".csv")
  write.csv(ech1, csv_name, row.names = FALSE, quote = F)
}  

###############################################################################################Ploidy data extraction and normalisation

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)

setwd("")

#
#####################
# controle haploide #
#####################

csv_ctrl2 = read.table('Control5.csv', header = T, sep = ',')
head(csv_ctrl2)
csv_ctrl2$ech <- as.character(csv_ctrl2$ech)
ctrl2_bin <- csv_ctrl2 %>% mutate(bin = GRN.B.HLog %/% 400)
sample = ctrl2_bin[1,1]

all_bins <- c(1:max(ctrl2_bin$bin)) %>% as.data.frame
colnames(all_bins) <- 'bin'

ctrl2_allbins <- full_join(ctrl2_bin, all_bins, by = 'bin') %>% mutate(ech = ifelse(is.na(ech), sample, ech)) %>% mutate(GRN.B.HLog = ifelse(is.na(ech), 0, GRN.B.HLog))

test <- ctrl2_allbins %>% group_by(ech, bin, GRN.B.HLog) %>% dplyr::summarise(count = length(GRN.B.HLog)) %>% as.data.frame
#test$GRN.B.HLog <- as.factor(test$GRN.B.HLog)
head(test)

#####################
# controle diploide #
#####################
csv_ctrl1 = read.table('Control8.csv', header = T, sep = ',')
csv_ctrl1$ech <- as.character(csv_ctrl1$ech)
ctrl1_bin <- csv_ctrl1 %>% mutate(bin = GRN.B.HLog %/% 400)
sample = ctrl1_bin[1,1]

ctrl1_allbins <- full_join(ctrl1_bin, all_bins, by = 'bin') %>% mutate(ech = ifelse(is.na(ech), sample, ech)) %>% mutate(GRN.B.HLog = ifelse(is.na(ech), 0, GRN.B.HLog))

test <- rbind(test,ctrl1_allbins %>% group_by(ech, bin, GRN.B.HLog) %>% dplyr::summarise(count = length(GRN.B.HLog)) %>% as.data.frame)
test_ctrl <- rbind(test,ctrl1_allbins %>% group_by(ech, bin, GRN.B.HLog) %>% dplyr::summarise(count = length(GRN.B.HLog)) %>% as.data.frame)
head(test)

#####################
# controle 3        #
#####################
csv_ctrl3 = read.table('Control3.csv', header = T, sep = ',')
csv_ctrl3$ech <- as.character(csv_ctrl3$ech)
ctrl3_bin <- csv_ctrl3 %>% mutate(bin = GRN.B.HLog %/% 400)
sample = ctrl3_bin[1,1]

ctrl3_allbins <- full_join(ctrl3_bin, all_bins, by = 'bin') %>% mutate(ech = ifelse(is.na(ech), sample, ech)) %>% mutate(GRN.B.HLog = ifelse(is.na(ech), 0, GRN.B.HLog))

test <- rbind(test,ctrl3_allbins %>% group_by(ech, bin, GRN.B.HLog) %>% dplyr::summarise(count = length(GRN.B.HLog)) %>% as.data.frame)
test_ctrl <- rbind(test,ctrl3_allbins %>% group_by(ech, bin, GRN.B.HLog) %>% dplyr::summarise(count = length(GRN.B.HLog)) %>% as.data.frame)
head(test)
tail (test)


###############
# Samples      #
###############

setwd("")
dossier = strsplit(getwd(), split = "/" )[[1]] %>% tail(n=1)

file_csv = dir(pattern='csv')

for(i in 1:length(file_csv)){
  #i=5
  print(i)
  
  csv_ech = read.table(file_csv[i], header = T, sep = ',')
  ech_bin <- csv_ech %>% mutate(bin = GRN.B.HLog %/% 400)
  file_name = strsplit(file_csv[i], '.', fixed = T)[[1]][1]
  strain_nb = strsplit(file_name,split = '_', fixed = T)[[1]] %>% tail(n=1)
  
  ech_allbins <- full_join(ech_bin, all_bins, by = 'bin') %>% mutate(ech = strain_nb) %>% mutate(GRN.B.HLog = ifelse(is.na(ech), 0, GRN.B.HLog))
  
  test <- rbind(test,ech_allbins %>% group_by(ech, bin, GRN.B.HLog) %>% dplyr::summarise(count = length(GRN.B.HLog)) %>% as.data.frame)
  
}

test$ech = as.factor(test$ech)
test = test %>% mutate(Hyb = ifelse(ech == "ctl5" | ech == "Haplo" | ech == "ctl3", "ctrl", "Ci"))
head(test)
tail(test)
nrow(test)


write.table(test, file="/tab_ploid_B3_wild.txt", quote = FALSE, sep = "\t", row.names = FALSE )

###################################################
Tab_P_heatmap = read.table("tab_ploid_B3_wild.txt",header=T,sep="\t")

head(Tab_P_heatmap)
Tab_P_heatmap = select (Tab_P_heatmap, -replicate)

#separate fluorescence values corresonding to 2n peaks from 4n peaks
test2n = filter (Tab_P_heatmap, GRN.B.HLog < 7500)
test4n = filter (Tab_P_heatmap, GRN.B.HLog > 7500)

#calculate maximum count cell for 2n fluorescence values
test_max2n = ddply(test2n, .(ech, Hyb), subset, count==max(count))

### Select samples different from 2n 
test_max2n0 = filter(test_max2n, count<20) # ploidy= 3n & 4n
pas2n = test_max2n0 %>% select(ech, Hyb)
testpas2n = left_join(pas2n, test4n)
test_max4n = ddply(testpas2n, .(ech, Hyb), subset,  count==max(count))

### Select diploid samples 2n 
test_max2ntrue= filter (test_max2n, count>20) # true diploids

###keep the 4n control
test_max4nC= ddply(test4n, .(ech, Hyb), subset, count==max(count))

max4nC = test_max4nC %>% select(ech, Hyb, count, GRN.B.HLog, bin)
test_max4nctl = filter(max4nC, ech== "Haplo", Hyb== "ctrl")

test_maxtot2n = full_join(test_max2ntrue, test_max4nctl)

#calculate the maximum cell count for each sample
GRN.B.HLog_max2n = test_maxtot2n %>% group_by(ech, Hyb) %>% dplyr::summarise(GRN.B.HLog_max2n = max(GRN.B.HLog)) %>% as.data.frame()
GRN.B.HLog_max4n = test_max4n %>% group_by(ech, Hyb) %>% dplyr::summarise(GRN.B.HLog_max4n = max(GRN.B.HLog)) %>% as.data.frame()

head (test_maxtot2n)
nrow(test_maxtot2n)

Tab_P_heatmap_max2n = left_join(GRN.B.HLog_max2n %>% mutate(GRN.B.HLog = GRN.B.HLog_max2n), test_maxtot2n) %>% unique() %>% mutate(bin_max = bin)
Tab_P_heatmap_max4n = left_join(GRN.B.HLog_max4n %>% mutate(GRN.B.HLog = GRN.B.HLog_max4n), test_max4n) %>% unique() %>% mutate(bin_max = bin)
Tab_P_heatmap_max2n = select (Tab_P_heatmap_max2n, -GRN.B.HLog_max2n)
Tab_P_heatmap_max4n = select (Tab_P_heatmap_max4n, -GRN.B.HLog_max4n)

Tab_P_heatmap_max = rbind (Tab_P_heatmap_max2n , Tab_P_heatmap_max4n)

head(Tab_P_heatmap_max)
Tab_P_heatmap$ech = as.character(Tab_P_heatmap$ech)
str(Tab_P_heatmap)

test_max1 = ddply(test, .(ech, Hyb), subset,  count==max(count))
test_max2 = select (test_max1, ech, bin, GRN.B.HLog, Hyb) 
head(test_max2, n = 100)
tail(test_max2)
summary(test_max2)
test_max2 %>% group_by(ech) %>% dplyr::summarise(nbligne=length(Hyb)) %>% as.data.frame()
test_max2 %>% filter(ech == 18) 

GRN.B.HLog_min = test_max2 %>% group_by(ech, Hyb) %>% dplyr::summarise(GRN.B.HLog_min = min(GRN.B.HLog)) %>% as.data.frame()
GRN.B.HLog_min %>% group_by(ech) %>% dplyr::summarise(nbligne=length(Hyb)) %>% as.data.frame()

#Normalisation 

GRN.B.HLog_min_normalised = mutate (Tab_P_heatmap_max, "normalisation"= 
                                                        ifelse (Hyb == "B", ((GRN.B.HLog / 3346.539)*1.122) - 0.117,
                                                              ifelse (Hyb == "C", ((GRN.B.HLog / 2864.384)*1.122) - 0.117,
                                                                      ifelse (Hyb == "A", ((GRN.B.HLog / 3409.676)*1.122) - 0.117, 
                                                                              ifelse (Hyb == "ctrl", ((GRN.B.HLog / 3766.7123)*1.122) - 0.117, "NA")))))

#save all tables for each cross type or wild strain separately
write.table(GRN.B.HLog_min_normalised, file="/tab_ploid_B3_wild.txt", quote = FALSE, sep = "\t", row.names = FALSE )



############################################################################################ Figure 1 B ploidy kinetics #######################################################

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)
library(lemon)
library(grid)
setwd("")

#read tables with fluorescence data for each cross and order them for the figure
data_P_all = read.table("Ploidy_all_4n_kinetics.txt",header=T,sep="\t")
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

pdf(("Ploidy_all_nbr_22_04_2020.pdf"), height = 10, width =9 )
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


############################################################################################ Figure 1 C ploidy frequency ########################################################

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)

setwd("")
Ploid_freq = read.table("Ploidy_frequency_total.txt",header=T,sep="\t")
blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

color<- c(red, red, red, "deepskyblue3", "deepskyblue3", "deepskyblue3", "maroon4", blue ,blue, green, green, black ,black)
Ploid_freq1 =filter(Ploid_freq, rep != "VL1")

pdf("ploidy_frequency_total_25_08_2020.pdf", width=16, height= 7)
ggplot(Ploid_freq, aes(x = (Hyb), y = as.numeric(WGD_freq), fill=Hyb)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=0.75)+ scale_x_discrete(limit = c("VL_S1", "VL_S2", "VL_A", "VL_C1","VL_C2", "VL_C3", "VL_B1","VL_B2", "VL_B3", "L1", "L2", "M1", "M2", "H1", "H2")) +
  scale_fill_manual(values=c("VL_B1" = red, "VL_B2" =red, "VL_B3" =red, "VL_C1" ="deepskyblue3", "VL_C2" ="deepskyblue3", "VL_C3" ="deepskyblue3", "VL_A" = "yellow4", "VL_S1" ="maroon4", "VL_S2" ="maroon4", "L1" =blue ,"L2" =blue, "M1" =green, "M2" =green, "H1" =black , "H2" =black))+
  ylab("Frequency of tetraploids") +
  xlab("Crosses")+
  theme(axis.text.x=element_text(size=18), axis.title.x=element_text(size=24)) +
  theme(axis.text.y=element_text(size=24), axis.title.y=element_text(size=24))
dev.off()


############################################################################# Figure S1 #########################################################################

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)
library(lemon)

setwd("")

#read tables with fluorescence data for each cross
table_Ploid_0 = read.table("tab_ploid_G_all.txt",header=T,sep="\t") %>% mutate(Hyb="VL_B3")
table_Ploid_1 = read.table("tab_ploid_I_all.txt",header=T,sep="\t") %>% mutate(Hyb="VL_B1")
table_Ploid_2 = read.table("tab_ploid_H_all.txt",header=T,sep="\t") %>% mutate(Hyb="VL_B2") 
table_Ploid_3 = read.table("tab_ploid_A_all.txt",header=T,sep="\t") %>% mutate(Hyb="L1") 
table_Ploid_4 = read.table("tab_ploid_D_all.txt",header=T,sep="\t") %>% mutate(Hyb="L2") 
table_Ploid_5 = read.table("tab_ploid_B_all.txt",header=T,sep="\t") %>% mutate(Hyb="M1") 
table_Ploid_6 = read.table("tab_ploid_E_all.txt",header=T,sep="\t") %>% mutate(Hyb="M2") 
table_Ploid_7 = read.table("tab_ploid_C_all.txt",header=T,sep="\t") %>% mutate(Hyb="H1") 
table_Ploid_8 = read.table("tab_ploid_F_all.txt",header=T,sep="\t") %>% mutate(Hyb="H2") 
table_Ploid_9 = read.table("tab_ploid_Hall_all.txt",header=T,sep="\t") %>% mutate(Hyb="VL_S2") 
table_Ploid_10 = read.table("tab_ploid_J_all.txt",header=T,sep="\t") %>% mutate(Hyb="VL_C1")
table_Ploid_11 = read.table("tab_ploid_K_all.txt",header=T,sep="\t") %>% mutate(Hyb="VL_C2")
table_Ploid_12 = read.table("tab_ploid_L_all.txt",header=T,sep="\t") %>% mutate(Hyb="VL_C3") 
table_Ploid_13 = read.table("tab_ploid_M_all.txt",header=T,sep="\t") %>% mutate(Hyb="VL_S1")
table_Ploid_14 = read.table("tab_ploid_N_all.txt",header=T,sep="\t") %>% mutate(Hyb="VL_A") 

#Remove empty wells
table_Ploid_0= filter(table_Ploid_0, ech!= 26)
table_Ploid_4= filter(table_Ploid_4, ech!= 93)

table_time = read.table("Info_Generations.txt",header=T,sep="\t")
table_time $Generations=as.factor(table_time $Generations)
tab_P_all= rbind (table_Ploid_0, table_Ploid_1, table_Ploid_2, table_Ploid_3, table_Ploid_4, table_Ploid_5, table_Ploid_6, table_Ploid_7, table_Ploid_8, table_Ploid_9, table_Ploid_10, table_Ploid_11, table_Ploid_12, table_Ploid_13, table_Ploid_14) 
tab_P_all= full_join(tab_P_all, table_time, by=c("Generations"))

tab_P_all1=filter(tab_P_all, normalisation<5)
tab_P_all2=filter(tab_P_all, normalisation>1)
tab_P_all= rbind(tab_P_all1, tab_P_all2)
tab_P_all= tab_P_all[!duplicated(tab_P_all), ]
write.table(tab_P_all, file="Tab_Ploidy_all_alltimepoint_06_2020.txt", quote = FALSE, sep = "\t", row.names = FALSE )

#keep only needed geerations 
tab_P_all=filter(tab_P_all, Generations== 22 | Generations== 352 | Generations==770 | Generations==1012)

write.table(tab_P_all, file="Tab_Ploidy_all_3timepoint_06_2020.txt", quote = FALSE, sep = "\t", row.names = FALSE )

tab_P_all = read.table("Tab_Ploidy_all_3timepoint_06_2020.txt",header=T,sep="\t")
tab_P_all = filter(tab_P_all, normalisation>1.4)

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


pdf(("Ploidy_3point_all_3n_12_06_2020.pdf"), height = 10, width =20)

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


############################################################################################### Figure S2 ######################################################################

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)
library(lemon)
library(grid)

setwd("")

#read tables with ploidy data for each cross
table_Ploid_1 = read.table("tab_ploid_JKL_4n_all.txt",header=T,sep="\t")
table_Ploid_2 = read.table("tab_ploid_BF4n_all2.txt",header=T,sep="\t") 
table_Ploid_3 = read.table("tab_ploid_AD4n_all.txt",header=T,sep="\t")
table_Ploid_41 = read.table("tab_ploid_ABDF_4n.txt",header=T,sep="\t")
Info_Ploid_4 = read.table("Info_ploid_ABDF_4n.txt",header=T,sep="\t")

Info_Ploid_4$ech=as.numeric(Info_Ploid_4$ech)
table_Ploid_41$ech=as.numeric(table_Ploid_41$ech)
table_Ploid_4 = full_join (table_Ploid_41, Info_Ploid_4, by=c("ech")) 

data_P_all_1= rbind (table_Ploid_1, table_Ploid_2, table_Ploid_3, table_Ploid_4) 

#Order the data for the figure
data_P_all_1 $Generations_f = factor(data_P_all_1 $Generations, levels=c(22, 88, 154, 220, 286, 352, 440, 460, 484, 550, 572, 616, 638, 682, 704 ,770))
data_P_all_1 $normalisation=as.numeric(data_P_all_1 $normalisation)
data_P_all_1= filter(data_P_all_1, Hyb_f != "NA")
data_P_all_1 $Generations=as.numeric(data_P_all_1 $Generations)

write.table(data_P_all_1, file="Tab_Ploidy_all_4n_05_2020.txt", quote = FALSE, sep = "\t", row.names = FALSE )

data_P_all_1 = read.table("Tab_Ploidy_all_4n_05_2020.txt",header=T,sep="\t")

data_P_all_1 $Hyb_f = factor(data_P_all_1 $Hyb, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))


blue_GC= rgb(0,0,(139/255))
red_GC= rgb((238/255),0,0)
green_GC= rgb(0,(139/255),0)
black_GC= rgb(0,0,0)


pdf(("Ploidy_4n_all_12_06_20.pdf"), height = 4, width =38 )
fig3B <- 
  ggplot(data_P_all_1, aes(x = Generations, y = as.numeric(normalisation), group=strain)) + 
  geom_line(alpha=2/10) + 
  geom_point(size=2.5, fill = "black", alpha=4/10) + 
  ylab("Ploidy")+ theme(axis.text=element_text(size=18))+
  xlab("")+
  facet_rep_grid(.~Hyb_f)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
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


################################################################################# Figure S4 ####################################################################
library("ggbeeswarm")

setwd("Ploidy-MA-CRL")

#read tables with fluorescence data for each cross
table_Ploid_1 = read.table("tab_ploid_B2_B3_wild_cor.txt",header=T,sep="\t")
table_Ploid_2 = read.table("tab_ploid_A_C_B1_wild_cor.txt",header=T,sep="\t") 

Tab_P_heatmap= rbind (table_Ploid_1, table_Ploid_2)

Tab_P_heatmap1= filter (Tab_P_heatmap, lineage == "SpB" | lineage == "SpC" | lineage == "SpA"| lineage == "Ctrl")
Tab_P_heatmap1$Lineage = factor(Tab_P_heatmap1$lineage, levels=c("Ctrl", "SpA", "SpB", "SpC"))

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)


pdf ('Wild_Ploidy_all_06_2020.pdf', width = 9 , height = 5)
ggplot(Tab_P_heatmap1, aes(x = Lineage, y = as.numeric(normalisation), fill = Lineage)) + geom_beeswarm(aes(color = Lineage, binwidth = 0.03)) +
  scale_color_manual(values=c("grey10" , "yellow4", red, "deepskyblue3"))+
  ylab("Ploidy")+ theme(axis.title =element_text(size=18))+
  xlab("Lineages")+ theme(axis.title=element_text(size=18))+
  theme(axis.text.x=element_text(size=18, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black"))+
  theme(axis.line = element_line(colour = "black"))
dev.off()


############################################################################################
#                                  Fertility data                                          #
############################################################################################ 

############################################################################################ Figure 2 B #################################################################

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
setwd("")
data = read.table("VL_tetraploid_fertility.txt", header = T, sep = '\t')
data $ploidy_f = factor(data $ploidy, levels=c("2n", "4n"))
data $cross_f = factor(data $cross, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))
data $time_f = factor(data $time, levels=c("Tini", "Tend"))


blue_GC= rgb(0,0,(139/255))
red_GC= rgb((238/255),0,0)
green_GC= rgb(0,(139/255),0)
black_GC= rgb(0,0,0)

data_4= filter(data, fertility != "NA")
data_4= filter(data_4, type == "tetraploid")


pdf(("VL_4n_fertility3_27_08_2020.pdf"), height = 5.5, width =13 )
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

############################################################################################ Figure S7 #################################################################
data_2= filter(data, fertility != "NA")
data_2= filter(data_2, type == "diploid")

pdf(("VL_4n_fertility_VL_2n_27_08_2020.pdf"), height = 3, width =8 )
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

library(tidyverse)
library(reshape2)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(magrittr)
library(stringr)
library(UpSetR)
library(gridExtra)
library(grid)
library(ggpubr)
library(gdata)
library(filesstrings)
library(cowplot)
library(gplots)

rm(list=ls())

############################################################### Calculate the maximum growth rate ###############################################################################

#read tables with OD data for each cross
setwd("")
data_all = read.table("Phenotypage_liq_raw_all_all_4n.txt", header = T, sep = '\t')

data_all = filter (data_all, strain != "NA")

#Transform  the time from second to hour
data_all = mutate (data_all, timeh =(time/3600))

#calculate slope

cal_slope2 <- function(size_vect, time_vect, window.regression = 10){
  #I create a temporary vector for the slopes and the intercepts
  slopes = NA
  max_slope_time = NA
  times = NA
  #it corresponds to the mean of the window
  fact_nas = sum(!is.na(size_vect))/length(size_vect)
  if (fact_nas > 0.8){
    for(j in 1:(length(size_vect) - window.regression + 1)){
      #I determine the start and end positions of the current window
      start = j
      end = j + window.regression - 1
      #I calculate the regression
      regression = lm(size_vect[start:end] ~ time_vect[start:end])
      slope = regression$coefficient[2]
      #I fill the temporary vectors
      slopes = c(slopes, slope)
      tim = as.numeric(time_vect[start])
      times = c(times, tim)
      #max_slope_time = c(max_slope_time, mean(time_vect[c(start,end)]))
    }
    df = data.frame("time"=times, "slope"=slopes)
    df = df[order(df$slope),]
    df = df[complete.cases(df), ]
    slope = df[length(df$time)-1,]
    #assemble regression parameters and choose the defined percentile for each

  }
  else {
    slope = NA

  }
  return(c(slope))

}

# Calculate the maximum slope
df_slopes2 = data_all %>% group_by(plate, sample, strain, cross, hyb, ploidy) %>%
  dplyr::summarise(max_slope = cal_slope2(OD, timeh)[[2]],
                   max_time = cal_slope2(OD, timeh)[[1]],
                   max_size = nth(OD, -3,order_by=OD),
                   ndiv = nth(OD, -3,order_by=OD)-min(OD, na.rm=T))

write.table(df_sum, file="Phenotypage_time_slope_2n.txt", sep = '\t')

data_all_slope = full_join(data_all, df_slopes2, by=c("plate", "sample", "strain", "cross", "hyb", "ploidy") )

#Remove max slopes that are after the second inflection of the growth curve due to the ADE2 deletion and absence of Adenine in the media
data_all_slope1 = filter(data_all_slope1, max_time < 50000)

data_F = filter(data_all_slope, (hyb=="A"| hyb=="L"))

pdf("curve_slope_time_point_L_2n.pdf", width=22, height= 18)

ggplot(data_F, aes(x = time, y = OD, group_by(starin, sample))) + 
  geom_smooth() + 
  geom_vline(data=data_F, aes(xintercept=max_time, colour="red"),linetype="dashed", size=1)+
  ylab("OD")+ theme(axis.text=element_text(size=11))+
  xlab("time")+
  facet_rep_wrap(.~strain)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'black', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
  theme(axis.text.x=element_text(size=16, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black")) #+
dev.off()
g <- ggplot_gtable(ggplot_build(fig3B))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c(blue_GC, blue_GC, blue_GC, green_GC, green_GC, 
           green_GC, green_GC, green_GC,blue_GC, blue_GC, black_GC, 
           black_GC, black_GC, black_GC, black_GC, red_GC, red_GC,
           red_GC, red_GC, red_GC, red_GC, red_GC, red_GC)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g) 
dev.off()


###########################################################################
data_all_slope1 = filter(data_all_slope, max_time < 48000)
df_sum <- data_all_slope1 %>% group_by(strain, cross, ploidy, hyb) %>%
  dplyr::summarize(median_max_slope = median(max_slope),
                   median_max_size = median(max_size),
                   sd_slope= sd(max_slope),
                   sd_size= sd(max_size))
write.table(df_sum, file="phenotyp_liq_growth_slope_max_time_all_4n_cor.txt", sep = '\t')


########################################################################################## Figure S8 #####################################################################


setwd("")
data1 = read.table("phenotyp_liq_growth_slope_max_time_all_4n_cor.txt", header = T, sep = '\t')
info1 = read.table("Info_phenotypage_all_4n.txt", header = T, sep = '\t')

data2 = read.table("phenotyp_liq_growth_slope_max_time_all_2n_2.txt", header = T, sep = '\t')
info2 = read.table("Info_phenotypage_all_2n.txt", header = T, sep = '\t')

data_all1 = full_join(data1, info1, by=c("strain"))%>% mutate("group"="tetraploid")

data = rbind(data2, data_all1)
data $time_f = factor(data $time, levels=c("Tini", "T2n", "T4n", "Tend"))
data=filter(data_all, time!="NA")
data $time_f = factor(data $time, levels=c("Tini", "T2n", "T4n", "Tend"))
data $hyb_f = factor(data $hyb, levels=c("VL3", "VL4", "VL5", "L1", "L2", "M1", "H2"))
data=filter(data, time_f!="NA")
data=filter(data, hyb_f!="NA")


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

pdf(("Phenotyp_liq_slope_max_all_all_2.pdf"), height = 4, width =20 )
fig3B <-
  ggplot(data, aes(x = time_f, y = as.numeric(median_max_slope), ymin=median_max_slope-sd_slope, ymax=median_max_slope+sd_slope, group=cross)) + 
  geom_line(aes(linetype= factor(group), alpha=2/10)) + 
  geom_pointrange(aes(colour = factor(group)), alpha=8/10) + 
  ylab("max_slope")+ theme(axis.text=element_text(size=18))+
  xlab("")+
  facet_rep_grid(.~hyb_f)+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+ theme(axis.title.y=element_text(size=18))+
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

############################################################ Calculate Fitness gain by WGD and post-WGD #############################################################################

library(ggplot2)
library(cowplot)
library(lattice)
library(gridExtra)
library(dplyr)
library(plyr)
library(lemon)
library(grid)
library(reshape2)
library(ggpubr)
library(formattable)
rm(list=ls())


setwd("")

data1 = read.table("phenotyp_liq_growth_slope_max_time_all_4n_cor.txt", header = T, sep = '\t')
info1 = read.table("Info_phenotypage_all_4n_2.txt", header = T, sep = '\t')
data2 = read.table("phenotyp_liq_growth_slope_max_time_all_2n_3.txt", header = T, sep = '\t')

data_all1 = full_join(data1, info1, by=c("strain"))%>% mutate("group"="tetraploid")
data_all1 =mutate(data_all1, "group1"="4n")
data2 =mutate(data2, "group1"="2n")
data = rbind(data2, data_all1)

#transform the maximum slopes by hour 
data =mutate(data, "max_slope_h"=median_max_slope*3600)
data =mutate(data, "sd_slope_h"=sd_slope*3600)

#Separate the data in groups before WGD and by WGD 
data_2n_4n=filter(data, group1 !="2n")
data_2n_4n=filter(data_2n_4n, time1=="T2n" |time1=="T4n")
data_2n_4n $time_f = factor(data_2n_4n $time1, levels=c("T2n", "T4n"))
data_2n_4n $hyb_f = factor(data_2n_4n $hyb, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))

data_2n_4n = filter(data_2n_4n , time_f !="NA" )
data_2n_4n = filter(data_2n_4n , hyb_f !="NA" )

#Separate the data post-WGD

data_4n_tend=filter(data, time1=="Tend" |time1=="T4n")
data_4n_tend=filter(data_4n_tend, group1 !="2n")

data_4n_tend=mutate(data_4n_tend, "slope_g"=(max_slope_h/Generations))
data_4n_tend=filter(data, time1=="Tend" |time1=="T4n" |time1=="Tmid")
data_4n_tend $time_f = factor(data_4n_tend $time1, levels=c("T4n", "Tend"))
data_4n_tend $hyb_f = factor(data_4n_tend $hyb, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))

compare_means(max_slope_h ~  time_f, method= "t.test",  data = data_2n_4n)

data $time_f = factor(data $time, levels=c("Tini", "T2n", "Tmid", "T4n", "Tend"))
data $hyb_f = factor(data $hyb, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))
data=filter(data, time_f!="NA")
data=filter(data, hyb_f!="NA")

data_T2n= filter(data, time=="T2n")
data_T4n= filter(data, time=="T4n" | time=="Tmid")
data_Tend= filter(data, time=="Tend")
data_Tini= filter(data, time=="Tini")

# Calculate the difference of max slope before and by WGD (T4n-T2n) and after WGD and at the end of the experiment  (Tend-T4n)
diff_2n_4n= left_join(data_T2n, data_T4n, by=c("cross", "hyb", "group1"))
diff_4n_tend= left_join(data_T4n, data_Tend, by=c("cross", "hyb", "group1"))
diff_tini_tend= left_join(data_Tini, data_Tend, by=c("cross", "hyb", "group1"))

diff_2n_4n=mutate(diff_2n_4n, "diff"= max_slope_h.y-max_slope_h.x)
diff_4n_tend=mutate(diff_4n_tend, "diff"= max_slope_h.y-max_slope_h.x)
diff_4n_tend=filter(diff_4n_tend, diff != "NA")
diff_tini_tend=mutate(diff_tini_tend, "diff"= max_slope_h.y-max_slope_h.x)
diff_tini_tend=filter(diff_tini_tend, diff != "NA")
diff_2n_4n=filter(diff_2n_4n, diff != "NA")

data_2=filter(data, group1 =="2n")
diff_2n_4n=filter(diff_2n_4n, group1 =="4n")

diff_4n_tend = mutate(diff_4n_tend, "nbr_Gener" =  (Generations.y-Generations.x))
diff_4n_tend = mutate(diff_4n_tend, "diff_Gener" =  (diff/nbr_Gener))

diff_4n_tend_4n=filter(diff_4n_tend, group1 =="4n")
diff_4n_tend_2n=filter(diff_4n_tend, group1 =="2n")

diff_tini_tend_4n= filter(diff_tini_tend, group1 =="4n")
diff_tini_tend_2n= filter(diff_tini_tend, group1 =="2n")

write.table (diff_4n_tend_4n, "tab_Fig_2D.txt", row.names = F, quote = F, sep= "\t")


########################################################################################### Figure S9 A ########################################################################

pdf(("Growth_rate_diff_2n_4n_27_08_2020.pdf"), height = 5, width =13 )
fig3B <-
  ggplot(data_2n_4n, aes(x = time_f, y = as.numeric(max_slope_h))) +#, group=hyb_f)) + 
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


########################################################################################### Figure S9 B ########################################################################

pdf("diff_4n_tend_h_4n_generations.pdf", width=14, height= 3.2)

fig3B <- ggplot(diff_4n_tend_4n[!is.na(diff_4n_tend$group1),], aes(x = nbr_Gener, y = as.numeric(diff), group=hyb)) +
  geom_point( shape=21, colour ="black", fill= "black", size=4, alpha= 0.6) +
  geom_hline(aes(yintercept = 0), colour = 'red', linetype = "dotted", size=1)+
  facet_grid(.~hyb_f.x, space= "free")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'white', face="bold", size = 18))+
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


########################################################################################### Figure 2 C ########################################################################

t_tests2 = diff_2n_4n %>%
  group_by(hyb_f.x) %>% dplyr::summarise(P = t.test(diff, mu = 0)$p.value, Sig = ifelse(P < 0.05, "*", "ns"), max_diff = max(diff))

t_tests2 <- t_tests2 %>%  mutate_if(is.numeric, round, digits = 2)

pdf("diff_2n_4n_h_normalised_27_08_2020.pdf", width=7.5, height= 5)
ggplot(diff_2n_4n, aes(x = hyb_f.x, y = as.numeric(diff), fill= hyb_f.x)) +
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


########################################################################################### Figure 2 D ########################################################################
diff_4n_tend_4n $hyb_f = factor(diff_4n_tend_4n $hyb, levels=c("VL_C1", "VL_C2", "VL_C3", "L1", "L2", "M1", "H2"))

diff_4n_tend_4n= filter (diff_4n_tend_4n, hyb != "L2" & hyb != "VL_C2") 

# Test whether each group differs from 0
t_tests = diff_4n_tend_4n %>% group_by(hyb) %>% dplyr::summarise(P = t.test(diff_Gener, mu = 0)$p.value, Sig = ifelse(P < 0.05, "*", "ns"), max_diff = max(diff_Gener))

t_tests <- t_tests %>%  mutate_if(is.numeric, round, digits = 3)



blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

pdf("diff_4n_tend_h_normalised_27_08_2020.pdf", width=7.5, height= 5)
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


########################################################################################### Calculate read coverage in windows of 10kb #####################

library(ggplot2)
library(dplyr)
library(reshape2)
rm(list=ls())

setwd("/MA_Sequences")
Table_info="Table_info_J.txt" 
Table_info <- read.table(Table_info, header = T, sep = "\t", fill = TRUE)

ref_cov_file="SpB_cov_J_strains.txt"  

ref_cov_wide <- read.table(ref_cov_file, header = T, sep = "\t", fill = TRUE)
ref_cov_long <- melt(ref_cov_wide, id.vars = c("contig", "pos"))
names(ref_cov_long)
head(ref_cov_long)

#Calculate the coverage over windows of 10 kb
depth_cov_10kb <- ref_cov_long %>% group_by(contig) %>% mutate(bin = pos %/% 10000) %>% group_by(contig, variable, bin) %>% summarise(cov = sum(value) / 10000) %>% mutate(Contig = gsub("_pilon", "", contig))  %>% as.data.frame()
depth_cov1_10kb =filter (depth_cov_10kb, bin != "NA")

ref = strsplit(ref_cov_file, "/", fixed = T)[[1]] %>% tail(1) %>% strsplit(., "_", fixed = T) %>% unlist() %>% head(1)

if (ref != "Scer") {
  ref_conversion_file = paste0(ref, "_tig_rearrangement.txt")
  
  ref_conversion <- read.table(ref_conversion_file, header = F, col.names = c("Chr", "contig", "Orientation"))
  
  depth_cov_10kb_conversion <- left_join(depth_cov1_10kb, ref_conversion, by = "contig")
  depth_cov_10kb_conversion$Chr_ord <- factor(depth_cov_10kb_conversion$Chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"))
  
} else {
  depth_cov_10kb_conversion <- depth_cov1_10kb %>% mutate(Chr_ord = gsub("_pilon", "", contig), Contig = Chr_ord)
}

depth_cov_10kb_conversion = filter(depth_cov_10kb_conversion, cov != "NA")

#Calculate the mean and median coverage by genome
summary_data <- depth_cov_10kb_conversion %>% group_by(variable) %>% summarise(mean_cov_tot = mean(cov), median_cov_tot = median(cov)) %>%  as.data.frame()

#Calculate the mean and median coverage by chromosome
summary_data_chr <- depth_cov_10kb_conversion %>% group_by(variable, contig) %>% summarise(mean_cov_chr = mean(cov), median_cov_chr = median(cov))  %>% as.data.frame()
tab_essai_meanChr = full_join(summary_data_chr, summary_data, by = "variable")
tab_essai_meanChr = full_join(tab_essai_meanChr, Table_info, by = "variable")
tab_essai = full_join(depth_cov_10kb_conversion, summary_data, by = "variable")
tab_essai = full_join(tab_essai, Table_info, by = "variable")


write.table (tab_essai, "/SpB_depth_cov_10kb_J_all.txt", row.names = F, quote = F, sep= "\t")
write.table (tab_essai_meanChr, "/SpB_meanChr_J_all.txt", row.names = F, quote = F, sep= "\t")


###################################################################### Figure S5 # and # Figure S12 #########################################################################################

library(ggplot2)
library(dplyr)
library(reshape2)
library(lattice)
library(gridExtra)
library(cowplot)
library(plyr)
library(lemon)
library(grid)
rm(list=ls())

setwd("")

#read table with coverages
d1 = read.table("SpB_depth_cov_10kb_Parents_all.txt", header=T) 
d2 = read.table("SpB_depth_cov_10kb_Parents_J.txt", header=T) 
d3 = read.table("SpB_depth_cov_10kb_Parents_K.txt", header=T) 
d4 = read.table("SpB_depth_cov_10kb_Parents_L.txt", header=T) 

Chr_corr = read.table("Chromosome_corr.txt", header=T) 

I = read.table("Parents_name.txt", header=T) 

#Join coverage and information tables
d1 = full_join(d1, Chr_corr, by=c("Chr", "Contig") )
d2 = full_join(d2, Chr_corr, by=c("Chr", "Contig") )
d3 = full_join(d3, Chr_corr, by=c("Chr", "Contig") )
d4 = full_join(d4, Chr_corr, by=c("Chr", "Contig") )


head(d12)
dall= rbind(d1, d2, d3, d4)
dall=select(dall, -Chr_ord)
#remove strains that are not considered in this study
dall=filter(dall, variable != "Jean.Talon")
dall=filter(dall, Chr != "NA")

write.table(dall, file="/Cov_all_Parents_10kb.txt", quote = FALSE, sep = "\t", row.names = FALSE )

dall$strain_ord <- factor(dall$variable, levels = c("MSH.604", "UWOPS.91.202", "LL2012_028", "LL2012_021", "LL2011_004", "LL2011_009", "MSH.587.1",  "LL2011_012", "LL2011_001", "YPS644", "YPS744", "LL2013_040", "LL2013_054"))

dall = filter(dall, strain_ord != "NA")
#Separate tables according to the timing
dall = filter(dall, strain_ord != "NA")
dall = filter(dall, Chr_cor != "12.5")
dall = filter(dall, Chr_cor != "15.5")

dall = full_join(dall, I, by=c("variable"))
dall$strain_ord <- factor(dall$parent, levels = c("SpB1", "SpB2", "SpB3", "SpB4", "SpC1", "SpC2", "SpC3",  "SpC4", "SpC5", "SpA1", "SpA2", "S.cer1", "S.cer2"))


############################################################################################################## Figure S12 ###############################################################
blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

pdf(file = paste0("Cov_Parents_all_lineblue_10kb.pdf"), height = 13, width = 10)

PlotA<-
  ggplot(dall, aes(x = bin, y = cov)) + geom_line() + geom_hline(data = dall, aes(yintercept = median_cov_tot), col = "blue3", size=1)+ 
  facet_grid(strain_ord~Chr_cor, scales = "free") + ylim(0, 325)+ 
  theme(axis.title.y = element_text(size=20), axis.text= element_text(), strip.text.y= element_text(face="bold", color= "white", size=18), strip.text.x= element_text(face="bold", color= "black", size=18))+
  labs(x=NULL) + labs(y='Read depth') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank())    

g <- ggplot_gtable(ggplot_build(PlotA))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c(red, red, red, red, "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "yellow4", "yellow4", "maroon4", "maroon4")

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

pdf(file = paste0("Cov_Parents_all_10kb.pdf"), height = 13, width = 10)

PlotA<-
  ggplot(dall, aes(x = bin, y = cov)) + geom_bar(stat="identity", position="identity", colour="grey30") + geom_hline(data = dall, aes(yintercept = median_cov_tot), col = "red", size=1)+ 
  facet_grid(strain_ord~Chr_cor, scales = "free") + ylim(0, 280)+ 
  theme(axis.title.y = element_text(size=20), axis.text= element_text(), strip.text.y= element_text(face="bold", color= "white", size=18), strip.text.x= element_text(face="bold", color= "black", size=18))+
  labs(x=NULL) + labs(y='Read depth') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank())    

g <- ggplot_gtable(ggplot_build(PlotA))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c(red, red, red, red, "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "yellow4", "yellow4", "maroon4", "maroon4")

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

############################################################################################################## Figure S5 ###############################################################

pdf(file = paste0("Cov_4n_tend_all_10kb_bar_cor.pdf"), height = 80, width = 40)

PlotA<-
  dall_Plot_tini <- ggplot(dall_Tend, aes(x = bin, y = cov)) + geom_bar(stat="identity", position="identity", colour="grey30") + geom_hline(data = dall_Tend, aes(yintercept = median_cov_tot), col = "red", size=2)+ 
  facet_grid(strain_ord~Chr_cor, scales = "free") + ylim(0, 280)+ 
  theme(axis.title.y = element_text(size=100), axis.text= element_text(), strip.text.y= element_text(face="bold", color= "white", size=40, margin = margin(0,1,0,1, "cm")), strip.text.x= element_text(face="bold", color= "black", size=80, margin = margin(1,0,1,0, "cm")))+
  labs(x=NULL) + labs(y='Read depth') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank())    

g <- ggplot_gtable(ggplot_build(PlotA))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3",
           "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", 
           "deepskyblue3", "deepskyblue3", "deepskyblue3", blue, blue, blue, blue, blue, green, green, green, green, black,
           black , black , black , black , black)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

pdf(file = paste0("Cov_4n_tend_all_10kb_bar.pdf"), height = 80, width = 40)

PlotA<-
  dall_Plot_tini <- ggplot(dall_Tend, aes(x = bin, y = cov)) + geom_bar() + geom_hline(data = dall_Tend, aes(yintercept = median_cov_tot), col = "red", size=2)+ 
  facet_grid(strain_ord~Chr, scales = "free") + ylim(0, 280)+ 
  theme(axis.title.y = element_text(size=100), axis.text= element_text(), strip.text.y= element_text(face="bold", color= "white", size=40, margin = margin(0,1,0,1, "cm")), strip.text.x= element_text(face="bold", color= "black", size=100, margin = margin(1,0,1,0, "cm")))+
  labs(x=NULL) + labs(y='Read depth') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank())    

g <- ggplot_gtable(ggplot_build(PlotA))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3",
           "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", 
           "deepskyblue3", "deepskyblue3", "deepskyblue3", blue, blue, blue, blue, blue, green, green, green, green, black,
           black , black , black , black , black)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

pdf(file = paste0("/Users/souhirmarsit/Desktop/MA_Sequences//","Cov_Parents_all_10kb_essai.pdf"), height = 55, width = 40)
PlotB<-
  dall_Plot_tend <- ggplot(dall_Tend, aes(x = pos, y = cov)) + geom_line() + geom_hline(data = dsall_Tend, aes(yintercept = median_cov), col = "red")+ 
  facet_grid(variable_f~chr_name, scales = "free") + ylim(0, 300)+ 
  theme(axis.title.y = element_text(size=70), axis.text= element_text(size=50), strip.text.y= element_text(face="bold", color= "white", size=73, margin = margin(0,1,0,1, "cm")), strip.text.x= element_text(face="bold", color= "black", size=80, margin = margin(1,0,1,0, "cm")))
+ labs(x=NULL) + labs(y='Read depth') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x=element_blank()) + theme(axis.ticks.x=element_blank())    

H <- ggplot_gtable(ggplot_build(PlotB))
stripr <- which(grepl('strip-r', H$layout$name))
fills <- c(blue, blue, blue, blue, blue, green, green, black)

l <- 1
for (z in stripr) {
  w <- which(grepl('rect', H$grobs[[z]]$grobs[[1]]$childrenOrder))
  H$grobs[[z]]$grobs[[1]]$children[[w]]$gp$fill <- fills[l]
  l <- l+1
}
grid.draw(H)
dev.off()



############################################################################################## Calculate Aneuploidy frequency ##############################################################

library(ggplot2)
library(dplyr)
library(reshape2)
rm(list=ls())


setwd("")
tab_essai_meanChr1 <- read.table("SpB_meanChr_F_4n.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr2 <- read.table("SpB_meanChr_B_4n.txt", header = T, sep = "\t", fill = TRUE)
tab_nomB <- read.table("nom_strain_B.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr2 = full_join(tab_essai_meanChr2, tab_nomB, by=c("variable"))
tab_essai_meanChr2 = filter(tab_essai_meanChr2, contig !="NA")
tab_essai_meanChr2 = filter(tab_essai_meanChr2, Hyb !="NA")
tab_essai_meanChr3 <- read.table("SpB_meanChr_J_4n.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr4 <- read.table("SpB_meanChr_K_4n.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr5 <- read.table("SpB_meanChr_L_4n.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr6 <- read.table("meanChr_F.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr7 <- read.table("meanChr_B.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr7 = filter(tab_essai_meanChr7, Hyb !="M1_40_Tend")
tab_essai_meanChr8 <- read.table("meanChr_A.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr9 <- read.table("meanChr_D.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr10 <- read.table("meanChr_C.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr11 <- read.table("meanChr_E.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr12 <- read.table("meanChr_H.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr13 <- read.table("meanChr_I.txt", header = T, sep = "\t", fill = TRUE)

tab_essai_meanChr = rbind(tab_essai_meanChr1, tab_essai_meanChr2, tab_essai_meanChr3, tab_essai_meanChr4, tab_essai_meanChr5, tab_essai_meanChr6, 
                          tab_essai_meanChr7, tab_essai_meanChr8, tab_essai_meanChr9, tab_essai_meanChr10, tab_essai_meanChr11, tab_essai_meanChr12, tab_essai_meanChr13)


Chr_corr <- read.table("chromosome_corr_freq.txt", header = T, sep = "\t", fill = TRUE)
tab_essai_meanChr = full_join(tab_essai_meanChr, Chr_corr, by=c("contig"))

#Calculate the copy number
tab_meanChr_norm = mutate (tab_essai_meanChr, "CNV"=((median_cov_chr/median_cov_tot)*Ploidy))

#Separate Tini and Tend data
cov_chr_Tini= filter (tab_meanChr_norm, time== "Tini")
cov_chr_Tend= filter (tab_meanChr_norm, time== "Tend")


################# Calculte the difference in copy number between Tend and Tini by Chromosome for tetraploids 4n
cov_chr_Tini= mutate (cov_chr_Tini, CNV=(CNV*2))

tab_cov_chr_all= full_join (cov_chr_Tini, cov_chr_Tend, by =c("strain_t", "cross", "contig", "Chr"))
tab_cov_chr_diff_A= mutate (tab_cov_chr_all, "CNV_diff"= ((CNV.y-CNV.x)))
tab_cov_chr_diff= mutate (tab_cov_chr_all, "CNV_diff"= ((CNV.y-CNV.x)))
tab_cov_chr_som0 = select (tab_cov_chr_diff_A, Chr_t.x, cross, Ploidy.x, Ploidy.y, CNV_diff , strain_t )
tab_cov_chr_som0= filter (tab_cov_chr_som0, CNV_diff != "NA")

################# Identify aneuploidies, chromosome gain and loss for  tetraploids 4n 
tab_cov_chr_som0_4n= filter (tab_cov_chr_som0, Ploidy.y == "4")
tab_cov_chr_som1_4n = mutate (tab_cov_chr_som0_4n, "Aneuploidy"=  
                                ifelse (CNV_diff > 0.75, 1,
                                        ifelse (CNV_diff < -0.75, 1, 0)))

tab_cov_chr_som1_gain_4n = mutate (tab_cov_chr_som0_4n, "Aneuploidy"=  
                                     ifelse (CNV_diff > 0.75, 1, 0))
tab_cov_chr_som1_loss_4n = mutate (tab_cov_chr_som0_4n, "Aneuploidy"=  
                                     ifelse (CNV_diff < -0.75, -1, 0))

aneup_gainloss = rbind(tab_cov_chr_som1_gain_4n, tab_cov_chr_som1_loss_4n)

write.table (aneup_gainloss, "Aneup_gain_loss_F_4n.txt", row.names = F, quote = F, sep= "\t")

################# Calculte the difference in copy number between Tend and Tini by Chromosome for diploids 2n

tab_cov_chr_all= full_join (cov_chr_Tini, cov_chr_Tend, by =c("strain_t", "cross", "contig", "Chr"))
tab_cov_chr_diff_A= mutate (tab_cov_chr_all, "CNV_diff"= ((CNV.y-CNV.x)))
tab_cov_chr_diff= mutate (tab_cov_chr_all, "CNV_diff"= ((CNV.y-CNV.x)))
tab_cov_chr_som0 = select (tab_cov_chr_diff_A, Chr_t.x, cross, Ploidy.x, Ploidy.y, CNV_diff , strain_t )
tab_cov_chr_som0= filter (tab_cov_chr_som0, CNV_diff != "NA")

################# Identify aneuploidies, chromosome gain and loss for diploids 2n 

tab_cov_chr_som0_2n= filter (tab_cov_chr_som0, Ploidy.y == "2")
tab_cov_chr_som1_2n = mutate (tab_cov_chr_som0_2n, "Aneuploidy"=  
                                ifelse (CNV_diff > 0.5, 1,
                                        ifelse (CNV_diff < -0.5, 1, 0)))

tab_cov_chr_som1_gain_2n = mutate (tab_cov_chr_som0_2n, "Aneuploidy"=  
                                     ifelse (CNV_diff > 0.5, 1, 0))
tab_cov_chr_som1_loss_2n = mutate (tab_cov_chr_som0_2n, "Aneuploidy"=  
                                     ifelse (CNV_diff < -0.5, -1, 0))

aneup_gainloss = rbind(tab_cov_chr_som1_gain_2n, tab_cov_chr_som1_loss_2n)
write.table (aneup_gainloss, "Aneup_gain_loss_J.txt", row.names = F, quote = F, sep= "\t")

################# Calculte the difference in copy number between Tend and Tini by Chromosome for triploids 3n

tab_cov_chr_all= full_join (cov_chr_Tini, cov_chr_Tend, by =c("strain_t", "cross", "contig", "Chr"))
tab_cov_chr_diff_A= mutate (tab_cov_chr_all, "CNV_diff"= ((CNV.y-CNV.x)))
tab_cov_chr_diff= mutate (tab_cov_chr_all, "CNV_diff"= ((CNV.y-CNV.x)))
tab_cov_chr_som0 = select (tab_cov_chr_diff_A, Chr_t.x, cross, Ploidy.x, Ploidy.y, CNV_diff , strain_t )
tab_cov_chr_som0= filter (tab_cov_chr_som0, CNV_diff != "NA")


################# Identify aneuploidies, chromosome gain and loss for triploids 3n

tab_cov_chr_som0_3n= filter (tab_cov_chr_som0, Ploidy.y == "3")
tab_cov_chr_som1_3n = mutate (tab_cov_chr_som0_3n, "Aneuploidy"=  
                                ifelse (CNV_diff > 0.6, 1,
                                        ifelse (CNV_diff < -0.6, 1, 0)))


tab_cov_chr_som1_gain_3n = mutate (tab_cov_chr_som0_3n, "Aneuploidy"=  
                                     ifelse (CNV_diff > 0.6, 1, 0))
tab_cov_chr_som1_loss_3n = mutate (tab_cov_chr_som0_3n, "Aneuploidy"=  
                                     ifelse (CNV_diff < -0.6, -1, 0))

aneup_gainloss = rbind(tab_cov_chr_som1_gain_3n, tab_cov_chr_som1_loss_3n)
write.table (aneup_gainloss, "Aneup_gain_loss_3n_A.txt", row.names = F, quote = F, sep= "\t")


######################################################################## Generate table of Aneuploidy frequency #############################################################

library(ggplot2)
library(dplyr)
library(reshape2)
library(grid)
rm(list=ls())

setwd("")
tab_essai_meanChr1 <- read.table("Aneup_gain_loss_A.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "2n") %>% mutate (Hybb = "L")
tab_essai_meanChr2 <- read.table("Aneup_gain_loss_3n_A.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "3n")%>% mutate (Hybb = "L")
tab_essai_meanChr3 <- read.table("Aneup_gain_loss_A_4n.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "4n")%>% mutate (Hybb = "L")
tab_essai_meanChr4 <- read.table("Aneup_gain_loss_D.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "2n")%>% mutate (Hybb = "L")
tab_essai_meanChr5 <- read.table("Aneup_gain_loss_3n_D.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "3n")%>% mutate (Hybb = "L")
tab_essai_meanChr6 <- read.table("Aneup_gain_loss_D_4n.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "4n")%>% mutate (Hybb = "L")
tab_essai_meanChr7 <- read.table("Aneup_gain_loss_B.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "2n")%>% mutate (Hybb = "M")
tab_essai_meanChr8 <- read.table("Aneup_gain_loss_B_4n.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "4n")%>% mutate (Hybb = "M")
tab_essai_meanChr9 <- read.table("Aneup_gain_loss_E.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "2n")%>% mutate (Hybb = "M")
tab_essai_meanChr10 <- read.table("Aneup_gain_loss_C.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "2n") %>% mutate (Hybb = "H")
tab_essai_meanChr11 <- read.table("Aneup_gain_loss_F.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "2n") %>% mutate (Hybb = "H")
tab_essai_meanChr12 <- read.table("Aneup_gain_loss_F_4n.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "4n")%>% mutate (Hybb = "H")
tab_essai_meanChr13 <- read.table("Aneup_gain_loss_H.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "2n") %>% mutate (cross = "VL_B2")%>% mutate (Hybb = "VL_B")
tab_essai_meanChr14 <- read.table("Aneup_gain_loss_I.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "2n") %>% mutate (cross = "VL_B1")%>% mutate (Hybb = "VL_B")
tab_essai_meanChr15 <- read.table("Aneup_gain_loss_J_4n.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "4n") %>% mutate (cross = "VL_C1")%>% mutate (Hybb = "VL_C")
tab_essai_meanChr16 <- read.table("Aneup_gain_loss_K_4n.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "4n") %>% mutate (cross = "VL_C2")%>% mutate (Hybb = "VL_C")
tab_essai_meanChr17 <- read.table("Aneup_gain_loss_L_4n.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "4n") %>% mutate (cross = "VL_C3")%>% mutate (Hybb = "VL_C")
tab_essai_meanChr18 <- read.table("Aneup_gain_loss_J.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "2n") %>% mutate (cross = "VL_C1")%>% mutate (Hybb = "VL_C")
tab_essai_meanChr19 <- read.table("Aneup_gain_loss_K.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "2n") %>% mutate (cross = "VL_C2")%>% mutate (Hybb = "VL_C")
tab_essai_meanChr20 <- read.table("Aneup_gain_loss_L.txt", header = T, sep = "\t", fill = TRUE) %>% mutate (ploidy = "2n") %>% mutate (cross = "VL_C3")%>% mutate (Hybb = "VL_C")


strain_tot_nmbr <- read.table("/strain_tot_nmbr.txt", header = T, sep = "\t", fill = TRUE)
chromID <- read.table("/Chromosome_id.txt", header = T, sep = "\t", fill = TRUE)


tab_tot_aneup = rbind(tab_essai_meanChr1, tab_essai_meanChr2, tab_essai_meanChr3, tab_essai_meanChr4, tab_essai_meanChr5, tab_essai_meanChr6, 
                      tab_essai_meanChr7, tab_essai_meanChr8, tab_essai_meanChr9, tab_essai_meanChr10, tab_essai_meanChr11, tab_essai_meanChr12, tab_essai_meanChr13, tab_essai_meanChr14,
                      tab_essai_meanChr15, tab_essai_meanChr16, tab_essai_meanChr17, tab_essai_meanChr18, tab_essai_meanChr19, tab_essai_meanChr20)

tab_tot_aneup1 = full_join(tab_tot_aneup, strain_tot_nmbr, by=c("cross", "ploidy"))
tab_tot_aneup = full_join(tab_tot_aneup1, chromID, by=c("Chr_t.x"))
tab_tot_aneup$total=as.numeric(tab_tot_aneup$total)
write.table (tab_tot_aneup, "/Aneuploidy_gainloss-by_chr_strain_all.txt", row.names = F, quote = F, sep= "\t")

tab_tot_aneup_gain = filter(tab_tot_aneup, Aneuploidy > 0| Aneuploidy == 0)%>% mutate (gainloss = "gain")
tab_tot_aneup_loss = filter(tab_tot_aneup, Aneuploidy < 0| Aneuploidy == 0)%>% mutate (gainloss = "loss")
tab_tot_aneup_gainloss = rbind (tab_tot_aneup_gain, tab_tot_aneup_loss)
tab_tot_aneup_gainloss1=unique(tab_tot_aneup_gainloss)
tab_tot_aneup_gainloss12 = filter (tab_tot_aneup_gainloss1, strain_t!= "D52" & strain_t!= "A74"& strain_t!= "D67")

write.table (tab_tot_aneup_gainloss12, "/Aneuploidy-by_chr_all_gainloss_25_01_2021.txt", row.names = F, quote = F, sep= "\t")

aneup_somme_gainloss1 = tab_tot_aneup_gainloss12  %>% group_by(strain_t, cross, ploidy, Hybb.x, Hybb.y, gainloss) %>% dplyr::summarise(Aneup_tot = sum(Aneuploidy)) %>% as.data.frame()
write.table (tab_tot_aneup, "/Aneuploidy-by_strain_all_gainloss.txt", row.names = F, quote = F, sep= "\t")


############################################################################################## Figure 3 A, Figure S10###########################

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


setwd("")

Aneup = read.table("Aneuploidy-by_strain_all_gainloss.txt", header = T, sep = "\t", fill = TRUE)


Aneup$cross_f <- factor(Aneup$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
Aneup$Hybb_f <- factor(Aneup$Hybb, levels = c("VL_C", "VL_B", "L","M", "H"))

Aneup= filter(Aneup, ploidy == "2n")
Aneup= filter(Aneup, Aneup_tot < 11)
Aneup= filter(Aneup, Aneup_tot > -11)
Aneup= filter(Aneup, Aneup_tot > -11)

Aneup_abs <- Aneup 
Aneup_abs$Aneup_tot <- abs(Aneup_abs$Aneup_tot)     

Aneup_abs_sum = Aneup_abs  %>% group_by(strain_t, cross, ploidy, Hyb, Hybb, cross_f, Hybb_f) %>% dplyr::summarise(Aneup_tot_tot = sum(Aneup_tot)) %>% as.data.frame()

Aneup_gain= filter(Aneup, gainloss == "gain")
Aneup_loss= filter(Aneup, gainloss == "loss")

Aneup_loss$Aneup_tot <- abs(Aneup_loss$Aneup_tot) 


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

Aneup.summary <- Aneup %>% group_by(cross_f, Hyb, Hybb_f, ploidy) %>% dplyr::summarise(mean=median( Aneup_tot), sd = sd(Aneup_tot, na.rm = TRUE)) 


pdf(file = paste0("Aneuploidy_rate_boxplot_corr_2n_25_08_2020_pvalue.pdf"), height = 6, width = 9)

my_comparisons <- list(c("VL_B", "VL_C"), c("VL_B", "L"), c("VL_B", "M"), c("VL_B", "H"),c("VL_C", "L"), c("VL_C", "M"), c("VL_C", "H"), c("L", "M"), c("L", "H"), c("M", "H"))
ggplot(Aneup_abs_sum[!is.na(Aneup_abs_sum$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y = 12.5, size= 4)+      # Add global p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)
   
my_comparisons <- list(c("VL_B", "VL_C"), c("VL_C", "M"), c("VL_C", "H"))
ggplot(Aneup_gain[!is.na(Aneup_gain$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y = 7.5, size= 4)+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)

my_comparisons <- list(c("VL_B", "VL_C"), c("VL_B", "L"), c("VL_B", "M"), c("VL_B", "H"), c("VL_C", "M"), c("VL_C", "H"), c("L", "H"))
ggplot(Aneup_loss[!is.na(Aneup_loss$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y =11, size= 4)+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)

dev.off()

############################################################################################## Figure S11 ##################################################
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


setwd("")

Aneup = read.table("Aneuploidy-by_strain_all_gainloss_Chr12.txt", header = T, sep = "\t", fill = TRUE)


Aneup$cross_f <- factor(Aneup$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
Aneup$Hybb_f <- factor(Aneup$Hybb, levels = c("VL_C", "VL_B", "L","M", "H"))

Aneup= filter(Aneup, ploidy == "2n")
Aneup= filter(Aneup, Aneup_tot < 11)
Aneup= filter(Aneup, Aneup_tot > -11)
Aneup= filter(Aneup, Aneup_tot > -11)

Aneup_abs <- Aneup 
Aneup_abs$Aneup_tot <- abs(Aneup_abs$Aneup_tot)     

Aneup_abs_sum = Aneup_abs  %>% group_by(strain_t, cross, ploidy, Hyb, Hybb, cross_f, Hybb_f) %>% dplyr::summarise(Aneup_tot_tot = sum(Aneup_tot)) %>% as.data.frame()

Aneup_gain= filter(Aneup, gainloss == "gain")
Aneup_loss= filter(Aneup, gainloss == "loss")

Aneup_loss$Aneup_tot <- abs(Aneup_loss$Aneup_tot) 


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

Aneup.summary <- Aneup %>% group_by(cross_f, Hyb, Hybb_f, ploidy) %>% dplyr::summarise(mean=median( Aneup_tot), sd = sd(Aneup_tot, na.rm = TRUE)) 


pdf(file = paste0("Aneuploidy_rate_boxplot_corr_2n_25_08_2020_Chr12.pdf"), height = 6, width = 9)
my_comparisons <- list(c("VL_B", "VL_C"), c("VL_B", "L"), c("VL_B", "M"), c("VL_B", "H"),c("VL_C", "L"), c("VL_C", "M"), c("VL_C", "H"), c("L", "M"), c("L", "H"), c("M", "H"))
ggplot(Aneup_abs_sum[!is.na(Aneup_abs_sum$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y = 8.5, size= 4)+      # Add global p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)

my_comparisons <- list(c("VL_B", "VL_C"), c("VL_C", "M"), c("VL_C", "H"), c("M", "H"))
ggplot(Aneup_gain[!is.na(Aneup_gain$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y = 7.5, size= 4)+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)

my_comparisons <- list(c("VL_B", "VL_C"), c("VL_B", "L"), c("VL_B", "M"), c("VL_B", "H"))
ggplot(Aneup_loss[!is.na(Aneup_loss$cross_f),], aes(x = Hybb_f, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3", red, red, blue, blue, green, green, black, black))+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(aes(group = Hybb_f),label.y =4, size= 4)+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)

dev.off()


############################################################################################### Figure 4 A, Figure S13 ####################################
rm(list=ls())

setwd("")

Aneup = read.table("Aneuploidy-by_strain_all_gainloss.txt", header = T, sep = "\t", fill = TRUE)

Aneup= filter(Aneup, cross != "VL_B1"& cross != "VL_B2"& cross != "M2"& cross != "H1")

Aneup$cross_f <- factor(Aneup$cross, levels = c("VL_C1", "VL_C2", "VL_C3", "VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
Aneup$Hybb_f <- factor(Aneup$Hybb, levels = c("VL_C", "VL_B", "L","M", "H"))

Aneup= filter(Aneup, ploidy == "2n"| ploidy == "4n")
Aneup= filter(Aneup, Aneup_tot < 11)
Aneup= filter(Aneup, Aneup_tot > -11)
Aneup= filter(Aneup, Aneup_tot > -11)

Aneup_abs <- Aneup 
Aneup_abs$Aneup_tot <- abs(Aneup_abs$Aneup_tot)     

Aneup_abs_sum = Aneup_abs  %>% group_by(strain_t, cross, ploidy, Hyb, Hybb, cross_f, Hybb_f) %>% dplyr::summarise(Aneup_tot_tot = sum(Aneup_tot)) %>% as.data.frame()

Aneup_gain= filter(Aneup, gainloss == "gain")
Aneup_loss= filter(Aneup, gainloss == "loss")

Aneup_loss$Aneup_tot <- abs(Aneup_loss$Aneup_tot) 

write.table (Aneup_abs_sum, "/tab_Fig_4A.txt", row.names = F, quote = F, sep= "\t")
write.table (Aneup_gain, "/tab_Fig_S13A.txt", row.names = F, quote = F, sep= "\t")
write.table (Aneup_loss, "/tab_Fig_S13B.txt", row.names = F, quote = F, sep= "\t")


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0,)

pdf(file = paste0("Aneuploidy_rate_boxplot_corr_4n_chr12_pvalue.pdf"), height = 6, width = 9)

my_comparisons <- list( c("2n", "4n"))

ggplot(Aneup_abs_sum[!is.na(Aneup_abs_sum$cross_f),], aes(x = ploidy, y = as.numeric(Aneup_tot_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3",  blue, blue, green, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4) 

ggplot(Aneup_abs_sum[!is.na(Aneup_abs_sum$cross_f),], aes(x = cross_f, y = as.numeric(Aneup_tot_tot), fill=ploidy)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3",  blue, blue, green, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means( method = "wilcox.test",  label = "p.format", size= 4)+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4) 

ggplot(Aneup_loss[!is.na(Aneup_loss$cross_f),], aes(x = ploidy, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3",  blue, blue, green, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4) 

ggplot(Aneup_loss[!is.na(Aneup_loss$cross_f),], aes(x = cross_f, y = as.numeric(Aneup_tot), fill=ploidy)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3",  blue, blue, green, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means( method = "wilcox.test",  label = "p.format", size= 4)

ggplot(Aneup_gain[!is.na(Aneup_gain$cross_f),], aes(x = ploidy, y = as.numeric(Aneup_tot), fill=cross_f)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3",  blue, blue, green, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4) 


ggplot(Aneup_gain[!is.na(Aneup_gain$cross_f),], aes(x = cross_f, y = as.numeric(Aneup_tot), fill=ploidy)) +
  geom_boxplot(alpha = 0.4)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.25)+
  scale_fill_manual(values=c("deepskyblue3", "deepskyblue3", "deepskyblue3",  blue, blue, green, black))+
  ylab("Aneuploidy / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means( method = "wilcox.test",  label = "p.format", size= 4)

dev.off()

###################################################################################Aneuploidy freq by chromosome heatmap###Figure 3 B, Figure 4 B #############

library(ggplot2)
library(dplyr)
library(reshape2)
library(grid)
library(cowplot)
library(lattice)
library(gridExtra)
library(plyr)
library(lemon)
#install.packages('wesanderson')
library(RColorBrewer)
library(wesanderson)
rm(list=ls())

tab_tot_aneup <- read.table("Aneuploidy-by_chr_all_gainloss_25_01_2021.txt", header = T, sep = "\t", fill = TRUE)

###### Calculate the frequency of chromosome gain or loss by cross by chromosome and by ploidy
tab_tot_aneup_2n = filter(tab_tot_aneup, ploidy=="2n")
tab_tot_aneup_3n = filter(tab_tot_aneup, cross=="L1" |cross=="L2")
tab_tot_aneup_4n = filter(tab_tot_aneup, ploidy=="4n")
tab_tot_aneup_2n$Hyb_f = factor(tab_tot_aneup_2n$Hybb.y, levels=c("H2_2n", "H1_2n", "M2_2n", "M1_2n", "L2_2n", "L1_2n", "VL_B2_2n", "VL_B1_2n", "VL_C3_2n", "VL_C2_2n", "VL_C1_2n"))
tab_tot_aneup_3n$Hyb_f = factor(tab_tot_aneup_3n$Hybb.y, levels=c(  "L1_2n", "L1_3n", "L1_4n", "L2_2n", "L2_3n", "L2_4n"))
tab_tot_aneup_4n$Hyb_f = factor(tab_tot_aneup_4n$Hybb.y, levels=c( "H2_4n", "M1_4n", "L2_4n", "L1_4n", "VL_C3_4n", "VL_C2_4n", "VL_C1_4n"))

tab_tot_aneup_2n_gain = filter(tab_tot_aneup_2n, gainloss=="gain")
tab_tot_aneup_3n_gain = filter(tab_tot_aneup_3n, gainloss=="gain")
tab_tot_aneup_4n_gain = filter(tab_tot_aneup_4n, gainloss=="gain")

tab_tot_aneup_2n_loss = filter(tab_tot_aneup_2n, gainloss=="loss")
tab_tot_aneup_3n_loss = filter(tab_tot_aneup_3n, gainloss=="loss")
tab_tot_aneup_4n_loss = filter(tab_tot_aneup_4n, gainloss=="loss")

tab_som_aneup_2n_gain = tab_tot_aneup_2n_gain %>% group_by(cross, ploidy, Chr, total_strain_nmbr, Hyb_f, gainloss) %>% dplyr::summarise(Aneup_freq_tot = sum(Aneuploidy)) %>% as.data.frame()
tab_som_aneup_3n_gain = tab_tot_aneup_3n_gain %>% group_by(cross, ploidy, Chr, total_strain_nmbr, Hyb_f, gainloss) %>% dplyr::summarise(Aneup_freq_tot = sum(Aneuploidy)) %>% as.data.frame()
tab_som_aneup_4n_gain = tab_tot_aneup_4n_gain %>% group_by(cross, ploidy, Chr, total_strain_nmbr, Hyb_f, gainloss) %>% dplyr::summarise(Aneup_freq_tot = sum(Aneuploidy)) %>% as.data.frame()

tab_som_aneup_2n_loss = tab_tot_aneup_2n_loss %>% group_by(cross, ploidy, Chr, total_strain_nmbr, Hyb_f, gainloss) %>% dplyr::summarise(Aneup_freq_tot = sum(Aneuploidy)) %>% as.data.frame()
tab_som_aneup_3n_loss = tab_tot_aneup_3n_loss %>% group_by(cross, ploidy, Chr, total_strain_nmbr, Hyb_f, gainloss) %>% dplyr::summarise(Aneup_freq_tot = sum(Aneuploidy)) %>% as.data.frame()
tab_som_aneup_4n_loss = tab_tot_aneup_4n_loss %>% group_by(cross, ploidy, Chr, total_strain_nmbr, Hyb_f, gainloss) %>% dplyr::summarise(Aneup_freq_tot = sum(Aneuploidy)) %>% as.data.frame()

tab_freq_aneup_2n_gain = mutate(tab_som_aneup_2n_gain, freq_aneup_chr=((Aneup_freq_tot/total_strain_nmbr)))
tab_freq_aneup_3n_gain = mutate(tab_som_aneup_3n_gain, freq_aneup_chr=((Aneup_freq_tot/total_strain_nmbr)))
tab_freq_aneup_4n_gain = mutate(tab_som_aneup_4n_gain, freq_aneup_chr=((Aneup_freq_tot/total_strain_nmbr)))

tab_freq_aneup_2n_loss = mutate(tab_som_aneup_2n_loss, freq_aneup_chr=((Aneup_freq_tot/total_strain_nmbr)))
tab_freq_aneup_3n_loss = mutate(tab_som_aneup_3n_loss, freq_aneup_chr=((Aneup_freq_tot/total_strain_nmbr)))
tab_freq_aneup_4n_loss = mutate(tab_som_aneup_4n_loss, freq_aneup_chr=((Aneup_freq_tot/total_strain_nmbr)))

tab_freq_aneup_2n_loss$freq_aneup_chr <- abs(tab_freq_aneup_2n_loss$freq_aneup_chr) 
tab_freq_aneup_3n_loss$freq_aneup_chr <- abs(tab_freq_aneup_3n_loss$freq_aneup_chr) 
tab_freq_aneup_4n_loss$freq_aneup_chr <- abs(tab_freq_aneup_4n_loss$freq_aneup_chr) 


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

######################################################################################## Figure 3 B ###################################################################################

pdf(file = paste0("Aneuploidy_freq_heatmap_legend_2n_25_08_2020.pdf"), height = 8, width = 20)

g_2n_g <- ggplot(tab_freq_aneup_2n_gain, aes(x = Chr , y = Hyb_f)) +
  geom_raster(aes (fill=freq_aneup_chr)) +  scale_x_discrete(expand = c(0, 0))+ scale_y_discrete(expand = c(0, 0))+
  scale_fill_gradient2(low = "white", high = "midnightblue", limits=c(0, 0.7)) + 
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25)) +
  ylab("") +xlab("")+ theme(legend.position = "bottom")
  coord_equal() 

g_2n_l <- ggplot(tab_freq_aneup_2n_loss, aes(x = Chr , y = Hyb_f)) +
  geom_raster(aes (fill=freq_aneup_chr)) +  scale_x_discrete(expand = c(0, 0))+ scale_y_discrete(expand = c(0, 0))+
  scale_fill_gradient2(low = "white", high = "darkred", limits=c(0, 0.7)) + 
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25)) +
  ylab("") +xlab("")+ theme(legend.position = "bottom")
  coord_equal() 
cowplot::plot_grid(g_2n_g, g_2n_l, ncol=2, rel_heights = c(300, 300), rel_widths = c(600, 600)) 

dev.off()

######################################################################################### Figure 4 B ##################################################################################

write.table (tab_freq_aneup_4n_gain, "/tab_Fig_4Bgain.txt", row.names = F, quote = F, sep= "\t")
write.table (tab_freq_aneup_4n_loss, "/tab_Fig_4Bloss.txt", row.names = F, quote = F, sep= "\t")


pdf(file = paste0("Aneuploidy_freq_heatmap_4n_25_08_2020.pdf"), height = 8, width = 20)
g_4n_g <- ggplot(tab_freq_aneup_4n_gain, aes(x = Chr , y = Hyb_f)) +
  geom_raster(aes (fill=freq_aneup_chr)) +  scale_x_discrete(expand = c(0, 0))+ scale_y_discrete(expand = c(0, 0))+
  scale_fill_gradient2(low = "white", high = "midnightblue", mid = "white", limits=c(0, 1)) + 
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25)) +
  ylab("") +xlab("")+
  coord_equal() 

g_4n_l <- ggplot(tab_freq_aneup_4n_loss, aes(x = Chr , y = Hyb_f)) +
  geom_raster(aes (fill=freq_aneup_chr)) +  scale_x_discrete(expand = c(0, 0))+ scale_y_discrete(expand = c(0, 0))+
   scale_fill_gradient2(low = "white", high = "darkred", limits=c(0, 1)) + 
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25)) +
  ylab("") +xlab("")+
  coord_equal()
cowplot::plot_grid(g_4n_g, g_4n_l, ncol=2, rel_heights = c(300, 300), rel_widths = c(600, 600)) 

dev.off()


################################################################################### Allele frequency table generation #########################################################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)
rm(list=ls())

setwd("")

#Read VCF merged table for each cross (here exemle for L cross)
d = read.table("var_L_snp.vcf_sep.txt", header=T, sep='\t') %>% mutate_all(as.character)
d = select(d, CHROM, POS, REF, ALT, contains("L"), contains("LL2011_009"), contains("LL2011_001"))
#Change the parents names according to each cross #MSH.604 #UWOPS.91.202 #LL2013_054 #LL2012_021

str(d)
names (d)
head (d)

#####filtrer for heterozygous variants between the two hybrid parents

d2 = mutate(d, "filtre_parent"=
              ifelse (LL2011_001.GT == LL2011_009.GT, "fail",
                      ifelse (LL2011_001.GT == "." | 	LL2011_009.GT == ".", "fail", "ok"))) 
head(d2)
str(d2)

d2$filtre_parent=as.factor(d2$filtre_parent)
summary(d2)

#filter d2 to keep only heterozygous variants

d3= filter (d2, filtre_parent == "ok")
head(d3)
summary(d3)

dtest= select (d3, CHROM , POS, REF, ALT, starts_with ("LL2011_001"), starts_with ("LL2011_009"), starts_with ("L"))

#Identify parental alleles

dtestlong = melt (dtest, id.vars = c("CHROM", "POS", "REF", "ALT"))

dtestvar= mutate(dtestlong, variable = gsub("AD.", "AD_", variable), strain = gsub("(^.*)[.].*$", "\\1", variable), var = gsub (".*[.](.*)", "\\1", variable))
head(dtestvar)
tail(dtestvar)

dtest=dcast(dtestvar, CHROM + POS + REF + ALT+ strain ~ var)

dtestallele= mutate(dtest, GT.AD_1 = gsub ("(.*)[/].*", "\\1", GT), GT.AD_2 = gsub (".*[/](.*)", "\\1", GT), GT.AD_0 = "B")
d_REF = dtestallele %>% filter(strain == "LL2011_001") %>%  select(-REF, -ALT, -starts_with("AD"), -GT, -GT.AD_2, -GT.AD_0) %>% dcast(CHROM + POS ~ strain)
head(d_REF)

dalBC= full_join(d_REF, dtestallele, by = c("CHROM", "POS"))

dalleleBC = mutate(dalBC, GT.AD_1= "C", GT.AD_2 = "MUT")

head(dalleleBC)

AD_0 = dalleleBC %>% select(CHROM, POS, strain, DP, AD, GT.AD_0) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_0")) %>% mutate(allele = GT.AD_0, AllelicDepth = value) %>% select(-GT.AD_0, -value)
head(AD_0)

AD_1 = dalleleBC %>% select(CHROM, POS, strain, DP, AD_1, GT.AD_1) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_1")) %>% mutate(allele = GT.AD_1, AllelicDepth = value) %>% select(-GT.AD_1, -value)
head(AD_1)

AD_2 = dalleleBC %>% select(CHROM, POS, strain, DP, AD_2, GT.AD_2) %>% melt(id.vars = c("CHROM", "POS", "strain", "DP", "GT.AD_2")) %>% mutate(allele = GT.AD_2, AllelicDepth = value) %>% select(-GT.AD_2, -value) %>% mutate(AllelicDepth = ifelse(AllelicDepth < 0, "NA", AllelicDepth), allele = ifelse( AllelicDepth == "NA", "NA", allele))
head(AD_2)

MONSUPERTABLEAU = rbind(AD_0, AD_1, AD_2)
tail(MONSUPERTABLEAU)

setwd("")

#Identification of contigs -> chromosome
ref_conversion <- read.table("SpC_tig_rearrangement.txt", header = F, col.names = c("Chr", "CHROM", "Orientation"))
head (ref_conversion)
SUPERTABLEAU <- left_join(MONSUPERTABLEAU, ref_conversion, by = "CHROM")
head (SUPERTABLEAU)
SUPERTABLEAU$Chr_ord <- factor(SUPERTABLEAU$Chr, levels = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"))

setwd("/Users/souhirmarsit/Desktop/MA_Sequences/All_genomes_LOH_analysis/")

Info_strain <- read.table("Table_info_L_LOH.txt", header = T, sep = "\t", fill = TRUE)
head (Info_strain)

SUPERTABLEAU$DP=as.integer(MONSUPERTABLEAU$DP)
SUPERTABLEAU$AllelicDepth=as.integer(MONSUPERTABLEAU$AllelicDepth)

#filter for variants with more than 20 reads 
tabDP10 = filter (SUPERTABLEAU, DP>20)

#Calculate allele frequency the ratio of allelic depth and the read depth of the locus   
tFA = mutate(tabDP10, "FA"= (AllelicDepth/DP))
head (tFA)

#filter for the allele frequency of Hybrid parent 2 (SpC or SpA or S. cerevisiae) of chr III
tFA_f = filter (tFA, allele == "C")
tFA_f <- left_join(tFA_f, Info_strain, by = "strain")
head(tFA_f)
tail(tFA_f)
summary(tFA_f)

write.table (tFA_f, "/tab_L_all_AF_C_05_2020.txt", row.names = F, quote = F, sep= "\t")


#################################################################################### LOH_density ### Figure S3 ##############################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)
library(cowplot)
library(lattice)
library(gridExtra)
library(lemon)
library(grid)
rm(list=ls())

setwd("")

I1 = read.table("Table_info_A_4n_AF.txt", header=T) 
I2 = read.table("Table_info_D_4n_AF.txt", header=T) 
I3 = read.table("Table_info_B_4n_AF.txt", header=T) 
I6 = read.table("Table_info_F_4n_AF.txt", header=T) 
I7 = read.table("Table_info_J_4n_AF.txt", header=T) 
I8 = read.table("Table_info_K_4n_AF.txt", header=T) 
I9 = read.table("Table_info_L_4n_AF.txt", header=T) 
I= rbind(I1, I2, I3, I6, I7, I8, I9)

d1 = read.table("tab_A_4n_AF.txt", header=T) 
d2 = read.table("tab_D_4n_AF.txt", header=T) 
d3 = read.table("tab_B40_49_4n_AF.txt", header=T) 
d4 = read.table("tab_F38_4n_AF.txt", header=T) 
d5 = read.table("tab_B_4n_AF_C_04_2020.txt", header=T) 
d6 = read.table("tab_F_4n_AF_C_04_2020.txt", header=T) 
d7 = read.table("tab_J_4n_AF.txt", header=T) 
d8 = read.table("tab_K_4n_AF.txt", header=T) 
d9 = read.table("tab_L_4n_AF.txt", header=T) 
d10 = read.table("tab_B_4n_AF.txt", header=T) 
d11 = read.table("tab_F_4n_AF.txt", header=T) 
d10= filter(d10, strain == "B32_P16"|strain == "B32_P35"| strain == "B59_P35")
d11= filter(d11, strain == "F47_P35")
# Remove strains that are repeated or not tetraploid
d5=filter(d5, strain != "B40_P1" & strain != "B40_P35" & strain != "B49_P1" & strain != "B49_P35" & strain != "B32_P16" & strain != "B32_P35" & strain != "B59_P35"& strain != "B55_P1" & strain != "B55_P35"  & strain != "B55_P16")
d6= filter(d6, strain != "F38_P1" & strain != "F38_P25" & strain != "F47_P35")

head(I)
head(d10)

d1 = select (d1, -time, -cross, -Ploidy, -Fertility, -Hyb)
d2 = select (d2, -time, -cross, -Ploidy, -Fertility)
d3 = select (d3, -time, -cross, -Ploidy, -Hyb)
d4 = select (d4, -time, -cross, -Ploidy, -Hyb)
d5 = select (d5, -time, -cross, -Ploidy, -Hyb)
d6 = select (d6, -time, -cross, -Ploidy)
d7 = select (d7, -time, -cross, -Ploidy, -Hyb)
d8 = select (d8, -time, -cross, -Ploidy, -Hyb)
d9 = select (d9, -time, -cross, -Ploidy, -Hyb)
d10 = select (d10, -time, -cross, -Ploidy, -Hyb)
d11 = select (d11, -time, -cross, -Ploidy, -Hyb)
d= rbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11)

head(d)
head(I)

dtot = full_join(d, I, by="strain")
head(dtot)
tail(dtot)

dtot1 = filter (dtot, time == "Tini"| time == "Tmid" | time == "Tend")
dtot1 = filter (dtot1, cross != "parent")
dtot2 = filter (dtot, cross == "parent")

write.table (dtot, "/tab_parents_4n_AF.txt", row.names = F, quote = F, sep= "\t")
dtot1 = read.table("/tab_parents_4n_AF.txt", header=T) 

#remove 3n strains
dtot1 = filter (dtot1, strain_t.y != "VLC1_40" & strain_t.y != "M1_55")
dtot1 = filter (dtot1, strain_t.y != "NA")

head(dtot1)
tail(dtot1)

dtot1$time_f = factor(dtot1$time, levels=c('Tini', 'Tmid', 'Tend'))
dtot1$strain_tt.y <- factor(dtot1$strain_t.y, levels = c("VLC1_2", "VLC1_13", "VLC1_17", "VLC1_27", "VLC1_32", "VLC1_48", "VLC1_53",  "VLC2_22", "VLC2_29", "VLC2_36", "VLC2_27", "VLC3_25", "VLC3_28", "VLC3_36" , "VLC3_43", "VLC3_50", "VLC3_57", "L1_31", "L1_51", "L1_87",
                                                         "L2_36", "L2_45", "M1_32", "M1_40", "M1_49",  "M1_55", "M1_59", "H2_18", "H2_38", "H2_43", "H2_47", "H2_57", "H2_61"))
#remove 2n B40 strain
dtot1 = filter(dtot1, strain != "B40_P35")

dtot_all = filter(dtot1, cross!= "L2")
dtot_L2 = filter(dtot1, cross== "L2")
dtot_L22 = filter(dtot_L2, FA > 0)
dtot1 = rbind (dtot_all, dtot_L22)

summary(dtot1)
str(dtot1)
dtot3= filter(dtot1, cross== "VL4")
dtot4= filter(dtot, strain== "LL2011_001" | strain== "LL2011_012"| strain== "MSH.587" | strain== "LL2011_004"| strain== "LL2011_009" | strain== "LL2013_054")
dtot5= mutate(dtot4, FA=(AllelicDepth/DP))
dtot3= filter(dtot, strain== "LL2013_054" )

blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

################################################################################# Figure S3 #####################################################################################

pdf(file = paste0("4n_all_density2_16_07_2020.pdf"), height = 20, width = 9)

PlotA<- ggplot(dtot1, aes(x = FA)) + geom_density()+ 
  facet_grid(strain_tt.y~ Ploidy, scales = "free") + 
  xlab("Allele frequency")+
  theme(axis.title.y = element_text(size=22), axis.title.x = element_text(size=22), axis.text.y= element_text(size=9), axis.text.x= element_text(size=20), strip.text.y= element_text(face="bold", color= "white", size=7.5))
  xlim(0, 1)

g <- ggplot_gtable(ggplot_build(PlotA))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3",
           "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3", "deepskyblue3",
           "deepskyblue3", "deepskyblue3", "deepskyblue3", blue, blue, blue, blue, blue, green, green, green, green, 
           black, black, black, black, black, black)

k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()


################################################################################### LOH_ heatmap ###### Fig S6 ###############################################################################

###########Analysis of VCF merged files of tetraploid strains 

library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)
library(lattice)
library(gridExtra)
library(cowplot)
library(lemon)

setwd("")


tFA_f = read.table("tab_All_4n_AF.txt", header = T, sep = "\t", fill = TRUE)
tFA_f_P = read.table("tab_F_P35_AF_C_04_2020.txt", header = T, sep = "\t", fill = TRUE)

tFA_f1 = filter (tFA_f, cross == "H2") 

tFA_f_P1 = filter (tFA_f_P, strain == "LL2013_054") 
tFA_f_P1 = mutate(tFA_f_P1 , FA_cor =(FA*1))

tFA_f_2n = filter (tFA_f1, Ploidy == 2 | Ploidy == 4) 
tFA_f_2n = mutate(tFA_f_2n , FA_cor =(FA*1))
tFA_f_3n = filter (tFA_f, Ploidy == 3)


######################################################
# Data of chromosome length
setwd("")

chr_length = read.table("chr_length_F.txt", header = T, sep = "\t", fill = TRUE)
chr_length= mutate(chr_length, "FA_cor"=1)

######################################################
setwd("All_genomes_LOH_analysis")
Info_strain <- read.table("Table_info_F_4n_AF.txt", header = T, sep = "\t", fill = TRUE)
head (Info_strain)

###################################################### Allele frquency for diploids 2n #####################################################################

tail(tFA_f_2n)
LOH_2n = select(tFA_f_2n, strain, CHROM, POS, FA_cor, allele) 
LOH_2n$POS = as.numeric(LOH_2n$POS)
#identify LOH, minimum of 3 successive SNPs with deviating allele frequency in a windows of 300bp
POS_100pb_2n <- LOH_2n %>% group_by(CHROM, strain) %>% mutate(bin = (POS %/% 300))%>% as.data.frame() 
POS_100pb2_2n <- POS_100pb_2n %>% add_count(CHROM, strain, bin, sort=TRUE)
POS_100pb3_2n = mutate(POS_100pb2_2n, "filtre"=
                         ifelse (FA_cor > 0.43  & FA_cor < 0.57 , "fail", "ok")) 
POS_100pb4_2n = filter(POS_100pb3_2n, n>3)
POS_100pb5_2n <- POS_100pb4_2n %>% add_count(CHROM, strain, bin, filtre, sort=TRUE) # %>% mutate(freq == nn/ sum(nn))
POS_100pb6_2n = mutate (POS_100pb5_2n, "SNP_freq_by_bin" = (nn/ n))
POS_100pb_LOH_2n = filter (POS_100pb6_2n, filtre == "ok")
POS_100pb_LOH1_2n = filter (POS_100pb_LOH_2n, SNP_freq_by_bin>0.85)

POS_100pb_LOH2_2n=full_join(chr_length, POS_100pb_LOH1_2n, by=c("CHROM", "POS", "strain", "FA_cor"))
setwd("/Users/souhirmarsit/Desktop/MA_Sequences")
ref_conversion <- read.table("SpB_tig_rearrangement.txt", header = F, col.names = c("Chr", "CHROM", "Orientation"))
head (ref_conversion)
POS_100pb_LOH2_2n <- left_join(POS_100pb_LOH2_2n, ref_conversion, by = "CHROM")
head (POS_100pb_LOH2_2n)
POS_100pb_LOH2_2n$Chr = as.numeric(POS_100pb_LOH2_2n$Chr)


LOH_f_2n <- full_join(POS_100pb_LOH2_2n, Info_strain, by = "strain")
tFA_f_P1 = select(tFA_f_P1, -DP, -variable, -AllelicDepth, -FA, -Chr_ord, -Fertility)
LOH_f_2n = select(LOH_f_2n, -bin, -n, -filtre, -nn, -SNP_freq_by_bin)
tFA_f_P1 = mutate(tFA_f_P1, "Hyb"="LL2013_054")
LOH_f_P = tFA_f_P1

LOH_f_Tini_2n = filter(LOH_f_2n, time == "Tini")
LOH_f_Tini_2n <- filter(LOH_f_Tini_2n, Ploidy != 3)
LOH_f_Tini_2n = filter (LOH_f_Tini_2n, Chr != "NA")
LOH_f_Tini_2n = rbind (LOH_f_P, LOH_f_Tini_2n)
LOH_f_Tini_2n$POS = as.numeric(LOH_f_Tini_2n$POS)
LOH_f_Tend_2n = filter(LOH_f_2n, time == "Tend"| time == "Tmid")
LOH_f_Tend_2n <- filter(LOH_f_Tend_2n, Ploidy != 3)
LOH_f_Tend_2n = rbind (LOH_f_P, LOH_f_Tend_2n)
LOH_f_Tend_2n$POS = as.numeric(LOH_f_Tend_2n$POS)
LOH_f_Tend_2n1 = rbind(LOH_f_Tend_2n, LOH_f_Tend_2n_47) 
LOH_f_Tend_47 = filter(LOH_f_2n, Hyb == "H2_47_P1")
LOH_f_Tend_end = filter (LOH_f_Tend_2n, Hyb == 'LL2013_054' | Hyb == 'H2_61_P35'| Hyb == 'H2_57_P35'| Hyb =='H2_43_P35'| Hyb =='H2_38_P16'| Hyb =='H2_18_P35')
LOH_f_Tend_mid = filter (LOH_f_Tend_2n, Hyb == 'LL2013_054' |Hyb == 'H2_61_P35'| Hyb == 'H2_57_P10'| Hyb =='H2_47_P35'| Hyb =='H2_43_P16'| Hyb =='H2_38_P16'| Hyb =='H2_18_P35')
LOH_f_Tend_end1 = rbind(LOH_f_Tend_end, LOH_f_Tend_47)

LOH_f_Tini_2n$Hyb_f = factor(LOH_f_Tini_2n$Hyb, levels=c('H2_18_P1', 'H2_38_P1', 'H2_43_P1', 'H2_47_P1', 'H2_57_P1', 'H2_61_P1', 'LL2013_054'))
LOH_f_Tend_end1$Hyb_f = factor(LOH_f_Tend_end1$Hyb, levels=c('H2_18_P35', 'H2_38_P16', 'H2_43_P35', 'H2_47_P1', 'H2_57_P35', 'H2_61_P35', 'LL2013_054'))
LOH_f_Tend_mid$Hyb_f = factor(LOH_f_Tend_mid$Hyb, levels=c('H2_18_P35', 'H2_38_P16', 'H2_43_P16', 'H2_47_P35', 'H2_57_P10', 'H2_61_P35', 'LL2013_054'))


pdf("rplot_AF_F_4n_Tend2", height = 4, width = 25)

g_tini <-
  ggplot(LOH_f_Tini_2n, aes(x= interaction(POS, Chr) , y = Hyb_f))+
  geom_raster(aes(fill=FA_cor)) + 
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1.08)) +
  facet_grid(~Chr, scales ="free_x", space="free_x")+
  theme_bw()+
  ylab("line")+ 
  theme(axis.text.y=element_text(size=18), axis.title.y=element_text(size=20)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text = element_text(colour = 'black', face="bold", size = 18))+ labs(x=NULL) + theme(axis.ticks=element_blank()) + theme(axis.text.x=element_blank())

g_tend <-
  ggplot(LOH_f_Tend_end1, aes(x= interaction(POS, Chr) , y = Hyb_f))+
  geom_raster(aes(fill=FA_cor)) + 
  scale_fill_gradient2(low = "darkred", high = "midnightblue", mid = "white", midpoint = 0.5, limits=c(0, 1.08)) +
  facet_grid(~Chr, scales ="free_x", space="free_x")+
  theme_bw()+
  ylab("line")+ 
  theme(axis.text.y=element_text(size=18), axis.title.y=element_text(size=18)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.text = element_text(colour = 'black', face="bold", size = 18))+ labs(x=NULL) + theme(axis.ticks=element_blank()) + theme(axis.text.x=element_blank())
+ gridExtra::grid.arrange(g_tini, g_tend, ncol=2)

dev.off()

#######################################################################################################identify the LOH_blocks#######################################

###########Analyse VCF merge toutes les souches
library(dplyr)
library(reshape2)
library(ggplot2)
library(plyr)
library(lattice)
library(gridExtra)
library(cowplot)
library(lemon)
rm(list=ls())

setwd("")

tFA_f = read.table("tab_A_AF_C_04_2020.txt", header = T, sep = "\t", fill = TRUE)

head(tFA_f)

tFA_f_P = filter (tFA_f, strain == "LL2011_004") 
tFA_f_P = mutate(tFA_f_P , FA_cor =(FA*1))
tFA_f_2n = filter (tFA_f, Ploidy == 2 ) 
tFA_f_2n = mutate(tFA_f_2n , FA_cor =(FA*1))
tFA_f_3n = filter (tFA_f, Ploidy == 3)


######################################################

setwd("")
chr_length = read.table("chr_length_A.txt", header = T, sep = "\t", fill = TRUE)

chr_length= mutate(chr_length, "FA_cor"=1)

######################################################
setwd("")
Info_strain <- read.table("Table_info_A_LOH.txt", header = T, sep = "\t", fill = TRUE)
head (Info_strain)

######################################################2n
tFA_f_2n= rbind (tFA_f_P, tFA_f_2n)
tail(tFA_f_2n)
LOH_2n = select(tFA_f_2n, strain, CHROM, POS, FA_cor, allele) 
LOH_2n$POS = as.numeric(LOH_2n$POS)
POS_100pb_2n <- LOH_2n %>% group_by(CHROM, strain) %>% mutate(bin = (POS %/% 300))%>% as.data.frame() 
POS_100pb2_2n <- POS_100pb_2n %>% add_count(CHROM, strain, bin, sort=TRUE)
POS_100pb3_2n = mutate(POS_100pb2_2n, "filtre"=
                         ifelse (FA_cor > 0.4  & FA_cor < 0.6 , "fail", "ok")) 
POS_100pb4_2n = filter(POS_100pb3_2n, n>3)
POS_100pb5_2n <- POS_100pb4_2n %>% add_count(CHROM, strain, bin, filtre, sort=TRUE) # %>% mutate(freq == nn/ sum(nn))
POS_100pb6_2n = mutate (POS_100pb5_2n, "SNP_freq_by_bin" = (nn/ n))
POS_100pb_LOH_2n = filter (POS_100pb6_2n, filtre == "ok")
POS_100pb_LOH1_2n = filter (POS_100pb_LOH_2n, SNP_freq_by_bin>0.85)
POS_100pb_LOH2_2n=full_join(chr_length, POS_100pb_LOH1_2n, by=c("CHROM", "POS", "strain", "FA_cor"))
setwd("/Users/souhirmarsit/Desktop/MA_Sequences")
ref_conversion <- read.table("SpB_tig_rearrangement.txt", header = F, col.names = c("Chr", "CHROM", "Orientation"))
head (ref_conversion)
POS_100pb_LOH2_2n <- left_join(POS_100pb_LOH2_2n, ref_conversion, by = "CHROM")
head (POS_100pb_LOH2_2n)
POS_100pb_LOH2_2n$Chr = as.numeric(POS_100pb_LOH2_2n$Chr)

LOH_f_2n <- left_join(POS_100pb_LOH2_2n, Info_strain, by = "strain")

POS_100pb_2n <- LOH_2n %>% group_by(CHROM, strain) %>% mutate(bin = (POS %/% 2000))%>% as.data.frame() 
POS_100pb2_2n <- POS_100pb_2n %>% add_count(CHROM, strain, bin, sort=TRUE)

LOH_f_P = filter (LOH_f_2n, strain=="LL2012_021") 
LOH_f_Tini_2n = filter(LOH_f_2n, time == "Tini")
LOH_f_Tini_2n <- filter(LOH_f_Tini_2n, Ploidy != 3)
LOH_f_Tini_2n = rbind (LOH_f_P, LOH_f_Tini_2n)
LOH_f_Tini_2n$POS = as.numeric(LOH_f_Tini_2n$POS)
LOH_f_Tend_2n = filter(LOH_f_2n, time == "Tend")
LOH_f_Tend_2n <- filter(LOH_f_Tend_2n, Ploidy != 3)
LOH_f_Tend_2n = rbind (LOH_f_P, LOH_f_Tend_2n)
LOH_f_Tend_2n$POS = as.numeric(LOH_f_Tend_2n$POS)
head(LOH_f_Tend_2n)
tail(LOH_f_Tend_2n)
LOH_f_Tend_2n_1 <- filter(LOH_f_Tend_2n,  FA_cor > 0.8)
LOH_f_Tend_2n_0 <- filter(LOH_f_Tend_2n,  FA_cor < 0.2)
LOH_f_Tend_2n_filter = rbind (LOH_f_Tend_2n_1, LOH_f_Tend_2n_0)
LOH_f_Tend_2n_filter1 <- LOH_f_Tend_2n_filter %>% add_count(CHROM, Chr, POS, sort=TRUE)

#Remove LOH that are present in more than 8 strains
LOH_f_Tend_2n_filter2 = filter(LOH_f_Tend_2n_filter1, nnn < 8)
LOH_f_Tend_2n_filter3 = filter(LOH_f_Tend_2n_filter2, strain != "F81_P35" & strain != "F38_P35" & strain != "F47_P35")

Aneup_by_chr = read.table("Aneup_by_chr_H_2n.txt", header = T, sep = "\t", fill = TRUE)

LOH_f_Tend_2n_filter3 = full_join(LOH_f_Tend_2n_filter3, Aneup_by_chr, by=c("strain_t", "cross", "Chr"))
LOH_f_Tend_2n_filter4 = filter(LOH_f_Tend_2n_filter3,  Aneuploidy != 1)
LOH_f_Tend_2n_filter5=unique(LOH_f_Tend_2n_filter4)

LOH_f_Tend_2n_filter5 <- LOH_f_Tend_2n_filter5 %>% add_count(CHROM, Chr, strain, bin, sort=TRUE)
LOH_f_Tend_2n_filter6 = filter(LOH_f_Tend_2n_filter5, nnnn > 5)

##############################################For aneuploid chromosomes

LOH_f_Tend_2n_filter1_aneup <- LOH_f_Tend_2n %>% add_count(CHROM, Chr, POS, sort=TRUE)
LOH_f_Tend_2n_filter2_aneup = filter(LOH_f_Tend_2n_filter1_aneup, nnn < 8)
LOH_f_Tend_2n_filter3_aneup = filter(LOH_f_Tend_2n_filter2_aneup, strain != "F81_P35" & strain != "F38_P35" & strain != "F47_P35")

Aneup_by_chr = read.table("Aneup_by_chr_H_2n.txt", header = T, sep = "\t", fill = TRUE)

LOH_f_Tend_2n_filter3_aneup = full_join(LOH_f_Tend_2n_filter3_aneup, Aneup_by_chr, by=c("strain_t", "cross", "Chr"))

LOH_f_Tend_2n_filter4_aneup = filter(LOH_f_Tend_2n_filter3_aneup,  Aneuploidy == 1)
LOH_f_Tend_2n_filter4_aneup1 <- LOH_f_Tend_2n_filter4_aneup %>% group_by(CHROM, Chr, strain) %>% dplyr::summarise(median_FA = median(FA_cor))%>% as.data.frame() 
LOH_f_Tend_2n_filter5_aneup = full_join(LOH_f_Tend_2n_filter4_aneup, LOH_f_Tend_2n_filter4_aneup1, by=c("strain", "CHROM", "Chr"))

LOH_f_Tend_2n_filter5_aneup1 = mutate(LOH_f_Tend_2n_filter5_aneup, "filtre_aneup"=
                                        ifelse ((FA_cor < (median_FA-0.15)), "ok", 
                                                ifelse ((FA_cor > (median_FA+0.15)), "ok", "fail"))) 

LOH_f_Tend_2n_filter5_aneup2 = filter (LOH_f_Tend_2n_filter5_aneup1, filtre_aneup == "ok")
LOH_f_Tend_2n_filter5_aneup3=unique(LOH_f_Tend_2n_filter5_aneup2)
LOH_f_Tend_2n_filter5_aneup3 <- LOH_f_Tend_2n_filter5_aneup3 %>% add_count(CHROM, Chr, strain, bin, sort=TRUE)
LOH_f_Tend_2n_filter6_aneup = filter(LOH_f_Tend_2n_filter5_aneup3, nnnn > 5)

#############################################################################join 2 tables

LOH_freq_1= select(LOH_f_Tend_2n_filter6, strain, Chr, POS, bin, FA_cor)
LOH_freq_aneup= select(LOH_f_Tend_2n_filter6_aneup, strain, Chr, POS, bin, FA_cor)

LOH_freq = rbind(LOH_freq_1, LOH_freq_aneup)

LOH_freq= filter (LOH_freq, strain!="NA")
LOH_freq= filter (LOH_freq, Chr!="NA")
LOH_freq= filter (LOH_freq, POS!="NA")
LOH_freq= filter (LOH_freq, bin!="NA")
LOH_freq= filter (LOH_freq, FA_cor!="NA")

LOH_freq1= filter(LOH_freq, FA_cor > 0)

LOH_f_P1 = select(LOH_f_P, strain, Chr, POS, bin, FA_cor)
LOH_freq_P = rbind (LOH_f_P1, LOH_freq1)


######################################################################################
write.table (LOH_freq, "/LOH_H_filter_aneup.txt", row.names = F, quote = F, sep= "\t")
###################################################################################
#Create blocks of LOH  
###############################

LOH_freq = read.table("/LOH_H_filter_aneup.txt", header = T, sep = "\t", fill = TRUE)

LOH_freq$bin=as.numeric(LOH_freq$bin)
LOH_freq$Chr=as.factor(LOH_freq$Chr)
LOH_freq$strain=as.factor(LOH_freq$strain)
LOH_freq$POS=as.numeric(LOH_freq$POS)
LOH_freq$FA_cor=as.numeric(LOH_freq$FA_cor)

LOH_freq$Chr <- as.factor(LOH_freq$Chr)

head(LOH_freq)
summary(LOH_freq)

# create an empty table that we fill using rbind
montableau <- NULL

# for each strain
for (Strain in levels(LOH_freq$strain)){
  print(paste("souche", Strain))
  
  db_strain <- filter(LOH_freq, strain == Strain)
 
  # Pour each chromosome
  for (chr in levels(db_strain$Chr)){
    print(paste("chr", chr))
    
    # extract a table db cotaining only the strain and the chromosome of interest 
    db <- filter(db_strain,  Chr == chr) %>% arrange_at("POS")
    
    if (nrow(db) > 0) {
      
      block <- 1
      blocklist <- 1
      
      # for each bin from db table
      for (i in c(2 : nrow(db))) {
  if ((db$bin[i]- db$bin[i-1]) <= 2){ block = block} else {block = block + 1}
        
       blocklist <- c(blocklist, block)
      }
      db$block <- blocklist
      
      montableau <- rbind(montableau, db)
      
    }
    
  }	
}

head(montableau, 100)
summary(montableau)
summary(LOH_freq)

montableau10 <- montableau %>% group_by(strain , Chr, block) %>% dplyr::summarise(median_FA = median(FA_cor)) %>% as.data.frame
montableau2 = select(montableau, -bin, -FA_cor)
montableau_all = full_join(montableau10, montableau2, by=c("strain", "Chr", "block"))

 #identif the parental origin
montableau1 = mutate(montableau_all, "origin"= ifelse ((median_FA < 0.22), "B",
                                                       ifelse ((median_FA > 0.8), "C", 
                                                               ifelse ((median_FA > 0.2 & median_FA < 0.4), "B2/3", 
                                                                       ifelse ((median_FA > 0.58 & median_FA < 0.75), "C2/3", "unknown")))))

LOH_block_start <- montableau1 %>% group_by(strain , Chr, block, origin) %>% dplyr::summarise(start = min(POS)) %>% as.data.frame
LOH_block_end <- montableau1 %>% group_by(strain , Chr, block, origin) %>% dplyr::summarise(end = max(POS)) %>% as.data.frame

LOH_block = full_join(LOH_block_start, LOH_block_end, by=c("strain", "Chr", "block", "origin"))
LOH_block = mutate(LOH_block, size=(end-start))
LOH_block_200pb = filter (LOH_block, size > 200)
LOH_block_500pb = filter (LOH_block, size > 500)
write.table (LOH_block_200pb, "/LOH_block_200pb_H_aneuploidy_cor.txt", row.names = F, quote = F, sep= "\t")
write.table (LOH_block_500pb, "/LOH_block_500pb_H_aneuploidy_cor.txt", row.names = F, quote = F, sep= "\t")
write.table (LOH_block_1000pb, "/LOH_block_1000pb_H_aneuploidy_cor.txt", row.names = F, quote = F, sep= "\t")



############################################################################################## Figure 3 C , D #############################################

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

setwd("")

LOH = read.table("LOH_segments_tot_08_2020.txt", header = T, sep = "\t", fill = TRUE)

LOH1=LOH

LOH1$cross_f <- factor(LOH1$cross, levels = c("VL_B1", "VL_B2", "L1", "L2","M1", "M2", "H1", "H2"))
LOH1$Hyb_f <- factor(LOH1$Hyb, levels = c("VL_B", "L","M", "H"))

LOH.summary <- LOH1 %>% group_by(cross_f, Hyb, ploidy) %>% dplyr::summarise(mean=mean(n), sd = sd(n, na.rm = TRUE)) 
LOH1= filter(LOH1, ploidy == "2n")
LOH1= filter(LOH1, Hyb == "L")
LOH1= filter(LOH1, ploidy == "2n"| ploidy == "3n")
LOH1= filter(LOH1, cross != "VL_B1"& cross != "VL_B2"&  cross != "M2"&  cross != "H1")


blue= rgb(0,0,(139/255))
red= rgb((238/255),0,0)
green= rgb(0,(139/255),0)
black= rgb(0,0,0)

pv <- LOH2 %>% group_by(Position, cross) %>% summarize(p.value = kruskal.test(n ~ cross_f)$p.value)

##################################################################################################### Figure 3 C #############################################


pdf(file = paste0("LOH_rate_boxplot_2n_Psignif_25_08_2020_cor.pdf"), height = 10, width = 10)

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


#################################################################################################### Figure 3 D #############################################

pdf(file = paste0("LOH_size_boxplot_2n_Psignif_08_2020.pdf"), height = 10, width = 10)

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

LOH1= filter(LOH1, Hyb == "L")
LOH1= filter(LOH1, ploidy == "2n"| ploidy == "3n")
LOH1= filter(LOH1, cross != "VL_B1"& cross != "VL_B2"&  cross != "M2"&  cross != "H1")
LOH2 <- LOH1 %>% add_count(strain, position, ploidy, sort=TRUE)


####################################################################################################### Figure 4 C ############################################

pdf(file = paste0("LOH_rate_boxplot_all_ploidy_point_4n2_08_2020_cor.pdf"), height = 7, width = 7)

my_comparisons <- list( c("2n", "4n"))
ggplot(LOH1[!is.na(LOH1$cross_f),], aes(x = ploidy, y = as.numeric(n), fill=ploidy)) +
  geom_boxplot(alpha = 0.4, color="black", size=0.25)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.35)+
  scale_fill_manual(values=c(blue, blue, blue, blue, green, green, black, black))+
  facet_grid(position~., space= "free")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'black', size = 20))+
  ylab("LOH / line")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)

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

pdf(file = paste0("LOH_size_boxplot_ploidy_4n2.pdf"), height = 7, width = 7)

my_comparisons <- list( c("2n", "4n"))
ggplot(LOH1[!is.na(LOH1$cross_f),], aes(x = ploidy, y = as.numeric(logsize), fill=ploidy)) +
  geom_boxplot(alpha = 0.4, color="black", size=0.25)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=0.35)+
  scale_fill_manual(values=c(blue, blue, green, black))+
  facet_grid(position~., space= "free")+ theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(strip.text = element_text(colour = 'black', size = 20))+
  ylab("log2(LOH size)")+ xlab("Cross")+ theme(axis.title.y = element_text(size=18), axis.title.x = element_text(size=18), axis.text.y= element_text(size=16), axis.text.x= element_text(size=18))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  label = "p.format", size= 4)


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


pdf(file = paste0("LOH_size_boxplot_ploidy_3n2_08_2020_cor.pdf"), height = 7, width = 7)

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

####################################
#               END                #
####################################




