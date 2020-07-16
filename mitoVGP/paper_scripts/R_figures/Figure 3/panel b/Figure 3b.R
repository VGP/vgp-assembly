setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(tibble)
library(ggpubr)
library(grid)
library(ggrepel)

table1 <- read.csv("heterozygosity_summary.txt", header=TRUE, sep="\t")

datalist = list()

for (i in 1:ncol(table1)) {
  data<-data.frame(cbind.data.frame(colnames(table1[i]), na.omit(table1[i])))
  colnames(data)<-c("group","value")
  datalist[[i]] <- data
}

HET = do.call(rbind, datalist)

t.test(HET$value, alternative = "greater")

table2 <- read.csv("regions.txt", header=FALSE, sep="\t")
colnames(table2)<-c("species","ID","start","end","type")
table2$match <- paste(table2$ID,"_",table2$type, sep="")
table2$diff <- (table2$end - table2$start)

sd<-table1 %>% lapply(na.omit) %>% lapply(sd)

sd<-do.call(rbind, lapply(sd, as.data.frame))
sd<-tibble::rownames_to_column(sd, "match")

data <- merge(sd,table2,by="match")

#CORRELATION repeat length/variance
cor_length_repeat <- data.frame(data$ID, data$diff, data$`X[[i]]`)
names(cor_length_repeat) <- c("dataset","repeat_size","sd")

cor(cor_length_repeat$repeat_size, cor_length_repeat$sd, method = "kendall")
cor(cor_length_repeat$repeat_size, cor_length_repeat$sd, method = "pearson")
cor(cor_length_repeat$repeat_size, cor_length_repeat$sd, method = "spearman")

ggscatter(cor_length_repeat, x = "repeat_size", y = "sd", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Repeat length (bp)", ylab = "Standard Deviation")

png("Fig. 3.png", width = 2000, height = 1200)

#plot

plot1<-ggplot(HET, aes(x = value, group = group, color = group)) +
  geom_area(aes(y = ..density.. * 53),
            size = 0.1, stat = 'density', alpha=0.01)+
  geom_point(data=HET,aes(y = group), position = position_jitter(w = 0, h = 0.2))+
  geom_vline(xintercept=0,
             color="black", linetype="dashed", size=1)+
  annotate("rect", xmin = -2250, xmax = -1350,  ymin = 0, ymax = Inf, alpha = .1)+
        scale_color_manual(values=rev(rainbow(ncol(table1), start = 0.1)))+
  annotate("rect", xmin = max(HET$value)+50, xmax = max(HET$value)+80,  ymin = 1, ymax = 17, alpha = .5)+
  annotate("rect", xmin = max(HET$value)+50, xmax = max(HET$value)+80,  ymin = 17+0.3, ymax = 23, alpha = .5)+
  annotate("rect", xmin = max(HET$value)+50, xmax = max(HET$value)+80,  ymin = 23+0.3, ymax = 27, alpha = .5)+
  annotate("rect", xmin = max(HET$value)+50, xmax = max(HET$value)+80,  ymin = 27+0.3, ymax = 46, alpha = .5)+
  annotate("rect", xmin = max(HET$value)+50, xmax = max(HET$value)+80,  ymin = 46+0.3, ymax = 68, alpha = .5)+
  annotate("text", x = max(HET$value)+60+70, y = (1+17)/2, label = "Fishes", family = "Arial",size = 15, hjust=0) +
  annotate("text", x = max(HET$value)+60+70, y = (17+23)/2+0.3, label = "Amphibians", family = "Arial",size = 15, hjust=0) +
  annotate("text", x = max(HET$value)+60+70, y = (23+27)/2+0.3, label = "Reptiles", family = "Arial",size = 15, hjust=0) +
  annotate("text", x = max(HET$value)+60+70, y = (27+46)/2+0.3, label = "Mammals", family = "Arial",size = 15, hjust=0) +
  annotate("text", x = max(HET$value)+60+70, y = (46+68)/2+0.3, label = "Birds", family = "Arial",size = 15, hjust=0)+
  theme(
    plot.margin = unit(c(1,15,1,1), "lines"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title = element_text(family = "Arial",
                              size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    #axis.ticks.y = element_blank(),
    axis.text.y=element_blank()
  )+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9), labels=function(x) format(x/1000, big.mark = ",", scientific = FALSE), limits = c(-2300,3900), expand=c(0,0))+
  labs(y="Repeat/Species", x = "Length deviation (kbp)")+
  geom_text_repel(aes(y=65,label=ifelse(value  == -1970,"C. anna","")), vjust=0, family = "Arial",size = 10, hjust=0, color="black", fontface = "italic", segment.color = NA)+
  geom_text_repel(aes(y=68,label=ifelse(value  == -2032,"H. comata","")), vjust=-2, hjust=-0.2, family = "Arial",size = 10, hjust=0, color="black", fontface = "italic", segment.color = NA)+
  geom_text_repel(aes(y=58,label=ifelse(value  == -2177,"C. canorus","")), vjust=0, family = "Arial",size = 10, hjust=0, color="black", fontface = "italic", segment.color = NA)+
  geom_text_repel(aes(y=21,label=ifelse(value  == -1464,"R. temporaria","")), vjust=1.6, hjust=1, family = "Arial",size = 10, hjust=0, color="black", fontface = "italic", segment.color = NA)+
  coord_cartesian(clip = 'off')
#  geom_text(aes(label=group, y=group),hjust=0, vjust=0)
  
plot1

dev.off()

