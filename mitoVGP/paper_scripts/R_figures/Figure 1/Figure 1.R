setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(dplyr)
library(scales)
library(tibble)
library(ggrepel)
library(ggpubr)
library(gtable)
library(grid)
library(readxl)

png("Fig. 1.png", width = 2000, height = 2000)

############NOVOPlasty###############

table1 <- read_excel("Summary - NOVOPlasty.xlsx", sheet = "Assembly_metadata")

#FIGURE 1a

is_outlier_id <- function(x) {
  return(x < 99.7)
}

table1<-table1 %>% filter(!table1$`VGP dataset` == "bCalAnn1")

out_identity_IUPAC <- data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`,identity = table1$`Identity (%)`, table1$Assembly)
out_identity_IUPAC <- out_identity_IUPAC %>% mutate(outlier = paste(sapply(table1$`Latin name`, substring, 1, 1), strsplit(as.character(table1$`Latin name`), " ") %>% sapply( "[", 2 ), sep=". ")) %>% mutate(is_outlier=ifelse(is_outlier_id(identity), identity, as.numeric(NA)))
out_identity_IUPAC$outlier[which(is.na(out_identity_IUPAC$is_outlier))] <- as.numeric(NA) 
out_identity_noIUPAC <- data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`,identity = table1$Rounding, table1$Assembly)

IDENTITY <- bind_rows("Identity" = out_identity_IUPAC, "IUPAC excluded" = out_identity_noIUPAC, .id = "group") %>% mutate(table1.Assembly = factor(table1.Assembly, levels=c("Multiple contigs", "Single contig", "Circular")))

#FIGURE 1b

#complete dataset
t.test(table1$`VGP length`,table1$`NOVOPlasty length`, paired = TRUE, alternative = "greater")
wilcox.test(table1$`VGP length`,table1$`NOVOPlasty length`, paired = TRUE, alternative = "greater")

table2 <- table1 %>% filter(Assembly == "Circular")

#Circular dataset
t.test(table2$`VGP length`,table2$`NOVOPlasty length`, paired = TRUE, alternative = "greater")
wilcox.test(table2$`VGP length`,table2$`NOVOPlasty length`, paired = TRUE, alternative = "greater")

is_outlier_ln <- function(x) {
  return(x < quantile(x, 0.20) - 1.5 * IQR(x) | x > quantile(x, 0.80) + 1.5 * IQR(x))
}

out_length <- data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`,length = table1$`NOVOPlasty length`, diff = table1$`VGP length` - table1$`NOVOPlasty length`, table1$Assembly)
out_length <- out_length %>% 
              mutate(outlier = paste(sapply(table1$`Latin name`, substring, 1, 1), strsplit(as.character(table1$`Latin name`), " ") %>%
              sapply( "[", 2 ), sep=". ")) %>% 
              mutate(is_outlier=ifelse(is_outlier_ln(diff), diff, as.numeric(NA))) %>%
              mutate(is_outlier_top10 = ifelse(is_outlier %in% tail(sort(is_outlier), 10), is_outlier, NA))
out_length$outlier[which(is.na(out_length$is_outlier_top10))] <- as.numeric(NA) 

VGP_length<-data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`, length = table1$`VGP length`, table1$Assembly)

LENGTH <- bind_rows("NOVOPlasty" = out_length, "VGP" = VGP_length, .id = "group") %>% mutate(table1.Assembly = factor(table1.Assembly, levels=c("Multiple contigs", "Single contig", "Circular")))

#FIGURE 1c

#complete dataset
t.test(table1$`VGP length`,table1$`NOVOPlasty length`, paired = TRUE, alternative = "greater")
wilcox.test(table1$`VGP length`,table1$`NOVOPlasty length`, paired = TRUE, alternative = "greater")

table3 <- table1 %>% filter(Assembly == "Circular")

#Circular dataset
t.test(table3$VGP_repeat,table3$NOVOplasty_repeat, paired = TRUE, alternative = "greater")
wilcox.test(table3$VGP_repeat,table3$NOVOplasty_repeat, paired = TRUE, alternative = "greater")

is_outlier_ln <- function(x) {
  return(x < quantile(x, 0.20) - 1.5 * IQR(x) | x > quantile(x, 0.80) + 1.5 * IQR(x))
}

out_repeat<- data.frame(table3$`VGP dataset`,table3$`Latin name`,table3$`Common name`,length = table3$VGP_repeat, diff = table3$VGP_repeat - table3$NOVOplasty_repeat, table3$Assembly)
out_repeat <- out_repeat %>%
              mutate(outlier = paste(sapply(table3$`Latin name`, substring, 1, 1), strsplit(as.character(table3$`Latin name`), " ") %>%
              sapply( "[", 2 ), sep=". ")) %>%
              mutate(is_outlier=ifelse(is_outlier_ln(diff), diff, as.numeric(NA))) %>%
              mutate(is_outlier_top10 = ifelse(is_outlier %in% tail(sort(is_outlier), 8), is_outlier, NA))
out_repeat$outlier[which(is.na(out_repeat$is_outlier_top10))] <- as.numeric(NA) 

NOVOPlasty_repeat<-data.frame(table3$`VGP dataset`,table3$`Latin name`,table3$`Common name`, length = table3$NOVOplasty_repeat, table3$Assembly)

REPEAT <- bind_rows("NOVOPlasty" = NOVOPlasty_repeat, "VGP" = out_repeat, .id = "group")

#MISSING

is_outlier_dm <- function(x) {
  return(x < quantile(x, 0.20, na.rm = TRUE) - 1.5 * IQR(x, na.rm = TRUE) | x > quantile(x, 0.80, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE))
}

out_missing <- data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`, missing = table1$`NOVOPlasty_missing count`, diff = table1$`NOVOPlasty_missing count` - table1$`VGP_missing count`)
out_missing <- out_missing %>% mutate(outlier = paste(sapply(table1$`Latin name`, substring, 1, 1), strsplit(as.character(table1$`Latin name`), " ") %>% sapply( "[", 2 ), sep=". ")) %>% mutate(is_outlier=ifelse(is_outlier_dm(diff), diff, as.numeric(NA)))
out_missing$outlier[which(is.na(out_missing$is_outlier))] <- as.numeric(NA) 

VGP_missing<-data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`, missing = table1$`VGP_missing count`)

MISSING <- bind_rows("NOVOPlasty" = out_missing, "VGP" = VGP_missing, .id = "group")

t.test(VGP_missing$missing,out_missing$missing,paired = TRUE, alternative = "less")
wilcox.test(VGP_missing$missing, out_missing$missing, paired = TRUE, alternative = "less")

###DUPS

out_duplicated <- data.frame(table3$`VGP dataset`,table3$`Latin name`,table3$`Common name`, duplicated = table3$`VGP_duplicated count`, diff = table3$`VGP_duplicated count` - table3$`NOVOPlasty_duplicated count`)
out_duplicated <- out_duplicated %>% mutate(outlier = paste(sapply(table3$`Latin name`, substring, 1, 1), strsplit(as.character(table3$`Latin name`), " ") %>% sapply( "[", 2 ), sep=". ")) %>% mutate(is_outlier=ifelse(is_outlier_dm(diff), diff, as.numeric(NA)))
out_duplicated$outlier[which(is.na(out_duplicated$is_outlier))] <- as.numeric(NA) 

NOVOPlasty_duplicated <- data.frame(table3$`VGP dataset`,table3$`Latin name`,table3$`Common name`, duplicated = table3$`NOVOPlasty_duplicated count`)

DUPLICATED <- bind_rows("NOVOPlasty" = NOVOPlasty_duplicated, "VGP" = out_duplicated, .id = "group")

t.test(out_duplicated$duplicated,NOVOPlasty_duplicated$duplicated,paired = TRUE, alternative = "greater")
wilcox.test(out_duplicated$duplicated, NOVOPlasty_duplicated$duplicated, paired = TRUE, alternative = "greater")

#plots

plot1<-ggplot(IDENTITY, aes(group, identity)) +
  geom_boxplot(width=0.5, size=1.5, aes(color=table1.Assembly), outlier.shape = NA, position = position_nudge(x = c(0))) +
  geom_point(aes(color=table1.Assembly), size=4) +
  geom_line(aes(group=table1..VGP.dataset.), colour="darkgray", size = 0.4, linetype = "dashed") +
  geom_label_repel(data = IDENTITY, aes(label = outlier), segment.size = 0.2, fontface = "italic",
                   nudge_x = -0.90, na.rm = TRUE, direction = "y", hjust = 0, label.size = 0, size=12, fill = NA)+
#  stat_compare_means(paired = TRUE, method.args = list(alternative = "greater"), label.x=1.5, size=12)+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title = element_text(family = "Arial",
                              size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    axis.title.x=element_blank(),
    plot.tag = element_text(size = rel(4), face = "bold")
  )+
  ylab("Identity (%)")+
  scale_color_manual(values=c("coral","gold","green3"))+
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  scale_x_discrete(expand = expansion(mult = c(0.95, 0.6))) +
  labs(tag = "a")

plot1

plot2<-ggplot(LENGTH, aes(group, length)) +
  geom_boxplot(width=0.5, size=1.5, aes(color=table1.Assembly), outlier.shape = NA, position = position_nudge(x = c(0))) +
  geom_point(aes(color=table1.Assembly), size=4, ) +
  geom_line(aes(group=table1..VGP.dataset.), colour="darkgray", size = 0.4, linetype = "dashed") +
  geom_label_repel(data = LENGTH, aes(label = outlier), segment.size = 0.2, fontface = "italic",
                   nudge_x = -1.45, na.rm = TRUE, direction = "y", hjust = 0, label.size = 0, size=12, fill = NA)+
  stat_compare_means(paired = TRUE, method.args = list(alternative = "greater"), label.x=0.75,label.y=25100, size=12)+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title = element_text(family = "Arial",
                              size = rel(4), colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    axis.title.x = element_blank(),
    plot.tag = element_text(size = rel(4), face = "bold")
  )+
  ylab("Assembly length (kbp)")+
  scale_color_manual(values=c("coral","gold","green3"))+
  scale_y_continuous(labels=function(x) format(x/1000, big.mark = ",", scientific = FALSE))+
  scale_x_discrete(expand = expansion(mult = c(1.30, 0.6))) +
  labs(tag = "b")

plot2

plot3<-ggplot(REPEAT, aes(group, length)) +
  geom_boxplot(width=0.5, size=1.5, aes(color=table3.Assembly), outlier.shape = NA, position = position_nudge(x = c(0))) +
  geom_point(aes(color=table3.Assembly), size=4, ) +
  geom_line(aes(group=table3..VGP.dataset.), colour="darkgray", size = 0.4, linetype = "dashed") +
  geom_label_repel(data = REPEAT, aes(label = outlier), segment.size = 0.2, fontface = "italic",
                   nudge_x = 0.90, na.rm = TRUE, direction = "y", hjust = 0, label.size = 0, size=12, fill = NA)+
  stat_compare_means(paired = TRUE, method.args = list(alternative = "greater"), label.x=1.25,label.y=3350, size=12)+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title = element_text(family = "Arial",
                              size = rel(4), colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    axis.title.x=element_blank(),
    plot.tag = element_text(size = rel(4), face = "bold")
  )+
  ylab("Repeat length (kbp)")+
  scale_color_manual(values=c("green3"))+
  scale_y_continuous(labels=function(x) format(x/1000, big.mark = ",", scientific = FALSE))+
  scale_x_discrete(expand = expansion(mult = c(0.6, 0.95))) +
  labs(tag = "c")

plot3

plot4<-ggplot(DUPLICATED, aes(group, duplicated)) +
  geom_boxplot(width=0.5, size=1.5, color="green3", outlier.shape = NA, position = position_nudge(x = c(0))) +
  geom_point(color="green3", size=4, ) +
  geom_line(aes(group=table3..VGP.dataset.), colour="darkgray", size = 0.4, linetype = "dashed") +
  geom_label_repel(data = DUPLICATED, aes(label = outlier), segment.size = 0.2, fontface = "italic",
                   nudge_x = 0.90, na.rm = TRUE, direction = "y", hjust = 0, label.size = 0, size=12, fill = NA)+
  stat_compare_means(paired = TRUE, method.args = list(alternative = "greater"), label.x=1.2,label.y=6, size=12)+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title = element_text(family = "Arial",
                              size = rel(4), colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    axis.title.x=element_blank(),
    plot.tag = element_text(size = rel(4), face = "bold")
  )+
  ylab("# gene duplications")+
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  scale_x_discrete(expand = expansion(mult = c(0.6, 0.95))) +
  labs(tag = "d")

plot4

##############################



###########Genbank############

is_outlier <- function(x) {
  return(x < quantile(x, 0.80) - 5 * IQR(x) | x > quantile(x, 0.20) + 5 * IQR(x))
}

table1 <- read_excel("Summary - Genbank.xlsx", sheet = "Assembly_metadata")

#FIGURE 1e

out_length <- data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`,length = table1$length1, diff = table1$length1 - table1$length2)
out_length <- out_length %>% mutate(outlier = paste(sapply(table1$`Latin name`, substring, 1, 1), strsplit(as.character(table1$`Latin name`), " ") %>% sapply( "[", 2 ), sep=". ")) %>% mutate(is_outlier=ifelse(is_outlier(diff), diff, as.numeric(NA)))
out_length$outlier[which(is.na(out_length$is_outlier))] <- as.numeric(NA) 

Genbank_length<-data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`, length = table1$length2)

LENGTH <- bind_rows("Genbank/RefSeq" = Genbank_length, "VGP" = out_length, .id = "group")

t.test(out_length$length,Genbank_length$length,paired = TRUE, alternative = "greater")
wilcox.test(out_length$length,Genbank_length$length, paired = TRUE, alternative = "greater")

####FIGURE 1f

out_repeat <- data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`,length = table1$VGP_repeat, diff = table1$VGP_repeat - table1$Genbank_repeat)
out_repeat <- out_repeat %>% mutate(outlier = paste(sapply(table1$`Latin name`, substring, 1, 1), strsplit(as.character(table1$`Latin name`), " ") %>% sapply( "[", 2 ), sep=". ")) %>% mutate(is_outlier=ifelse(is_outlier(diff), diff, as.numeric(NA)))
out_repeat$outlier[which(is.na(out_repeat$is_outlier))] <- as.numeric(NA) 

Genbank_repeat<-data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`, length = table1$Genbank_repeat)

REPEAT <- bind_rows("Genbank/RefSeq" = Genbank_repeat, "VGP" = out_repeat, .id = "group")

t.test(out_repeat$length,Genbank_repeat$length,paired = TRUE, alternative = "greater")
wilcox.test(out_repeat$length,Genbank_repeat$length, paired = TRUE, alternative = "greater")

####FIGURE 1g
out_missing <- data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`, missing = table1$`Genbank_missing count`, diff = table1$`Genbank_missing count` - table1$`VGP_missing count`)
out_missing <- out_missing %>% mutate(outlier = paste(sapply(table1$`Latin name`, substring, 1, 1), strsplit(as.character(table1$`Latin name`), " ") %>% sapply( "[", 2 ), sep=". ")) %>% mutate(is_outlier=ifelse(is_outlier(diff), diff, as.numeric(NA)))
out_missing$outlier[which(is.na(out_missing$is_outlier))] <- as.numeric(NA) 

VGP_missing<-data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`, missing = table1$`VGP_missing count`)

MISSING <- bind_rows("Genbank/RefSeq" = out_missing, "VGP" = VGP_missing, .id = "group")

t.test(VGP_missing$missing,out_missing$missing,paired = TRUE, alternative = "less")
wilcox.test(VGP_missing$missing, out_missing$missing, paired = TRUE, alternative = "less")

####FIGURE 1h
out_duplicated <- data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`, duplicated = table1$`VGP_duplicated count`, diff = table1$`VGP_duplicated count` - table1$`Genbank_duplicated count`)
out_duplicated <- out_duplicated %>% mutate(outlier = paste(sapply(table1$`Latin name`, substring, 1, 1), strsplit(as.character(table1$`Latin name`), " ") %>% sapply( "[", 2 ), sep=". ")) %>% mutate(is_outlier=ifelse(is_outlier(diff), diff, as.numeric(NA)))
out_duplicated$outlier[which(is.na(out_duplicated$is_outlier))] <- as.numeric(NA) 

Genbank_duplicated <- data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`, duplicated = table1$`Genbank_duplicated count`)

DUPLICATED <- bind_rows("Genbank/RefSeq" = Genbank_duplicated, "VGP" = out_duplicated, .id = "group")

t.test(out_duplicated$duplicated,Genbank_duplicated$duplicated,paired = TRUE, alternative = "greater")
wilcox.test(out_duplicated$duplicated, Genbank_duplicated$duplicated, paired = TRUE, alternative = "greater")

#plots

plot5<-ggplot(LENGTH, aes(group, length)) +
  geom_boxplot(width=0.5, size=1.5, aes(color=group), outlier.shape = NA, position = position_nudge(x = c(0))) +
  geom_point(aes(color=group), size=4) +
  geom_line(aes(group=table1..VGP.dataset.), colour="darkgray", size = 0.4, linetype = "dashed") +
  geom_label_repel(data = LENGTH, aes(label = outlier), segment.size = 0.2, fontface = "italic",
                   nudge_x = 0.90, na.rm = TRUE, direction = "y", hjust = 0, label.size = 0, size=12, fill = NA)+
  stat_compare_means(paired = TRUE, method.args = list(alternative = "greater"), label.x=1.25, label.y=25000, size=12)+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title = element_text(family = "Arial",
                              size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    axis.title.x=element_blank(),
    plot.tag = element_text(size = rel(4), face = "bold")
  )+
  ylab("Assembly length (kbp)")+
  scale_color_manual(values=c("#335eaa","#ffd700"))+
  scale_y_continuous(labels=function(x) format(x/1000, big.mark = ",", scientific = FALSE))+
  scale_x_discrete(expand = expansion(mult = c(0.6, 0.95))) +
  labs(tag = "e")


plot5

plot6<-ggplot(REPEAT, aes(group, length)) +
  geom_boxplot(width=0.5, size=1.5, aes(color=group), outlier.shape = NA, position = position_nudge(x = c(0))) +
  geom_point(aes(color=group), size=4) +
  geom_line(aes(group=table1..VGP.dataset.), colour="darkgray", size = 0.4, linetype = "dashed") +
  geom_label_repel(data = REPEAT, aes(label = outlier), segment.size = 0.2, fontface = "italic",
                   nudge_x = 0.90, na.rm = TRUE, direction = "y", hjust = 0, label.size = 0, size=12, fill = NA)+
  stat_compare_means(paired = TRUE, method.args = list(alternative = "greater"), label.x=1.25, label.y=5000, size=12)+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title = element_text(family = "Arial",
                              size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    axis.title.x=element_blank(),
    plot.tag = element_text(size = rel(4), face = "bold")
  )+
  ylab("Repeat length (kbp)")+
  scale_color_manual(values=c("#335eaa","#ffd700"))+
  scale_y_continuous(labels=function(x) format(x/1000, big.mark = ",", scientific = FALSE))+
  scale_x_discrete(expand = expansion(mult = c(0.6, 0.95))) +
  labs(tag = "f")

plot6

plot7<-ggplot(MISSING, aes(group, missing)) +
  geom_boxplot(width=0.5, size=1.5, aes(color=group), outlier.shape = NA, position = position_nudge(x = c(0))) +
  geom_point(aes(color=group), size=4) +
  geom_line(aes(group=table1..VGP.dataset.), colour="darkgray", size = 0.4, linetype = "dashed") +
  geom_label_repel(data = MISSING, aes(label = outlier), segment.size = 0.2, fontface = "italic",
                   nudge_x = -0.90, na.rm = TRUE, direction = "y", hjust = 0, label.size = 0, size=12, fill = NA)+
  stat_compare_means(paired = TRUE, method.args = list(alternative = "less"), label.x=1.25, label.y=5, size=12)+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title = element_text(family = "Arial",
                              size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    axis.title.x=element_blank(),
    plot.tag = element_text(size = rel(4), face = "bold")
  )+
  ylab("# missing genes")+
  scale_color_manual(values=c("#335eaa","#ffd700"))+
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  scale_x_discrete(expand = expansion(mult = c(0.6, 0.95))) +
  labs(tag = "g")

plot7

plot8<-ggplot(DUPLICATED, aes(group, duplicated)) +
  geom_boxplot(width=0.5, size=1.5, aes(color=group), outlier.shape = NA, position = position_nudge(x = c(0))) +
  geom_point(aes(color=group), size=4) +
  geom_line(aes(group=table1..VGP.dataset.), colour="darkgray", size = 0.4, linetype = "dashed") +
  geom_label_repel(data = DUPLICATED, aes(label = outlier), segment.size = 0.2, fontface = "italic",
                   nudge_x = 0.90, na.rm = TRUE, direction = "y", hjust = 0, label.size = 0, size=12, fill = NA)+
  stat_compare_means(paired = TRUE, method.args = list(alternative = "greater"), label.x=1.25, label.y=6, size=12)+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title = element_text(family = "Arial",
                              size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    axis.title.x=element_blank(),
    plot.tag = element_text(size = rel(4), face = "bold")
  )+
  ylab("# gene duplications")+
  scale_color_manual(values=c("#335eaa","#ffd700"))+
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  scale_x_discrete(expand = expansion(mult = c(0.6, 0.95))) +
  labs(tag = "h")

plot8

g1 <- ggplotGrob(plot1)
g2 <- ggplotGrob(plot2)
g3 <- ggplotGrob(plot3)
g4 <- ggplotGrob(plot4)
g5 <- ggplotGrob(plot5)
g6 <- ggplotGrob(plot6)
g7 <- ggplotGrob(plot7)
g8 <- ggplotGrob(plot8)

gr1 <- cbind(g1, g2, size = "first")
gr2 <- cbind(g3, g4, size = "first")
gr3 <- cbind(g5, g6, size = "first")
gr4 <- cbind(g7, g8, size = "first")

g <- rbind(gr1,gr2,gr3,gr4)
gtable_add_padding(g, unit(1, "cm"))
colnames(g) <- letters[1:8]
grid.newpage()
grid.draw(g)

dev.off()

####SUPP FIGURES

#CORRELATION length/repeats

cor_length_repeat <- data.frame(out_length$table1..VGP.dataset., table1$Status, out_length$diff, out_repeat$diff)
names(cor_length_repeat) <- c("dataset","status","lengths_diff","repeats_diff")

cor(cor_length_repeat$lengths_diff, cor_length_repeat$repeats_diff, method = "kendall")
cor(cor_length_repeat$lengths_diff, cor_length_repeat$repeats_diff, method = "pearson")
cor(cor_length_repeat$lengths_diff, cor_length_repeat$repeats_diff, method = "spearman")

png("../Extended Data Fig. 5/Extended Data Fig. 5.png", width = 2000, height = 2000)

ggscatter(cor_length_repeat, x = "repeats_diff", y = "lengths_diff",
          size = 10, cor.coef.size = 20,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Repeat length difference (bp)", ylab = "Assembly length difference (bp)")+
  theme(
    plot.margin = margin(t = 20, b = 0),
    axis.title = element_text(family = "Arial",
                              size = rel(6),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(6), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
  )+
  scale_x_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))

dev.off()

##GC

out_gc <- data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`, gc = table1$gc1, diff = table1$gc1 - table1$gc2)
out_gc <- out_gc %>% mutate(outlier = paste(sapply(table1$`Latin name`, substring, 1, 1), strsplit(as.character(table1$`Latin name`), " ") %>% sapply( "[", 2 ), sep=". ")) %>% mutate(is_outlier=ifelse(is_outlier(diff), diff, as.numeric(NA)))
out_gc$outlier[which(is.na(out_gc$is_outlier))] <- as.numeric(NA) 

Genbank_gc<-data.frame(table1$`VGP dataset`,table1$`Latin name`,table1$`Common name`, gc = table1$gc2)

GC <- bind_rows("Genbank/RefSeq" = Genbank_gc, "VGP" = out_gc, .id = "group")

png("../Extended Data Fig. 6/Extended Data Fig. 6.png", width = 2000, height = 700)

ggplot(GC, aes(group, gc)) +
  geom_boxplot(width=0.5, size=1.5, aes(color=group), outlier.shape = NA) +
  geom_point(aes(color=group), size=4) +
  geom_line(aes(group=table1..VGP.dataset.), colour="darkgray", size = 0.4, linetype = "dashed") +
  #geom_label_repel(data = GC, aes(label = outlier), segment.size = 0.2, fontface = "italic",
  #                 nudge_x = 0.90, na.rm = TRUE, direction = "y", hjust = 0, label.size = 0, size=12, fill = NA)+
  stat_compare_means(paired = TRUE, method.args = list(alternative = "two.sided"), label.x=1.25, label.y=50, size=12)+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.position = "none",
    axis.title = element_text(family = "Arial",
                              size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    axis.title.x=element_blank(),
    plot.tag = element_text(size = rel(4), face = "bold")
  )+
  ylab("GC content (%)")+
  scale_color_manual(values=c("#335eaa","#ffd700"))+
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  scale_x_discrete(expand = expansion(mult = c(0.6, 0.95)))

dev.off()


