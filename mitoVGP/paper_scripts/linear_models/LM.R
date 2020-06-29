# load libraries
library(dplyr)
library(ggplot2)
library(lawstat)
library(readxl)
library(lmPerm)
library(multcomp)

# set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# import data

mitoVGP.All <- read_excel("../Supplementary Tables/Supplementary Tables.xlsx", sheet = "ST1")

names(mitoVGP.All)
# select relevant variables for analyses
mitoVGP <- mitoVGP.All[,c(1, 7, 8, 9, 11, 13, 15, 17, 25, 28)]

# set factors and numbers
mitoVGP <- mitoVGP %>% mutate_if(is.character, as.factor)

mitoVGP <- mutate(mitoVGP, Size.selection.cutoff = ifelse(Size.selection.kbp.Pacbio > 20,
                                                          "Above 20", "Below 20"))
str(mitoVGP)
# mitoVGP$Taxonomic.group <- as.factor(mitoVGP$Taxonomic.group)
# mitoVGP$Tissue.type.Pacbio <- as.factor(mitoVGP$Tissue.type.Pacbio)
# mitoVGP$DNA.extraction.Pacbio <- as.factor(mitoVGP$DNA.extraction.Pacbio)
# mitoVGP$Library.prep.fragmentation.Pacbio <- as.factor(mitoVGP$Library.prep.fragmentation.Pacbio)
# mitoVGP$Library.prep.Pacbio <- as.factor(mitoVGP$Library.prep.Pacbio)


# produce and save result tables
out1 <- with(mitoVGP, table(Tissue.type.Pacbio, Taxonomic.group))
out2 <- with(mitoVGP, table(Tissue.type.Pacbio, Success))
out3 <- with(mitoVGP, table(Size.selection.kbp.Pacbio, Success))
out4 <- with(mitoVGP, table(Size.selection.cutoff, Success))

write.table(out1, file='../Supplementary Tables/ST3.tsv', quote=FALSE, sep='\t')
write.table(out2, file='../Supplementary Tables/ST5.tsv', quote=FALSE, sep='\t')
write.table(out3, file='../Supplementary Tables/ST6a.tsv', quote=FALSE, sep='\t')
write.table(out4, file='../Supplementary Tables/ST6b.tsv', quote=FALSE, sep='\t')

# code for Figure 3 in the extended data

mitoVGP <- mitoVGP %>%
  mutate(Success = ifelse(Success == "no",0,1))

png("../Extended Data Fig. 3/Extended Data Fig. 3.png", width = 2000, height = 750)

ggplot(mitoVGP, aes(x=reorder(Tissue.type.Pacbio, -Available.Pacbio.mtDNA.reads),  
                    y=Available.Pacbio.mtDNA.reads)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=8,
               outlier.size=4)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("")+
  ylab("Number of Pacbio reads")+
  theme(
    plot.margin = margin(t = 20, b = 0),
    axis.title = element_text(family = "Arial",
                              size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6))
  )

dev.off()

# other plots
plot(Success ~ Available.Pacbio.mtDNA.reads, data = mitoVGP)
boxplot(Available.Pacbio.mtDNA.reads ~ Taxonomic.group, data = mitoVGP)

# mod.mitoVGP <- glm(Success ~ `Available Pacbio mtDNA reads`, family="binomial", data=mitoVGP)
# summary(mod.mitoVGP)
# anova(mod.mitoVGP, test="Chi")


png("../Extended Data Fig. 2/Extended Data Fig. 2.png", width = 2000, height = 750)

mitoVGP %>%
  mutate(size = ifelse(Available.Pacbio.mtDNA.reads < 1, "0", ">0")) %>%
  ggplot(aes(x = Available.Pacbio.mtDNA.reads, y = Success, colour = size))  +
  stat_sum() +
  geom_smooth(method = "glm", formula = y ~ x, se = FALSE,
              method.args = list(family = binomial()), colour = "gray")+
  scale_x_continuous(trans='sqrt', breaks = c(5,100,200,300,400,500))+
  scale_color_manual(values = c('green', 'red'))+
  scale_y_continuous(name="Success", breaks = c(0, 1), labels = c("Failure", "Success"))+
  scale_size_continuous(range = c(5, 15))+
  xlab("Number of Pacbio reads")+
  theme(
    plot.margin = margin(t = 20, b = 0),
    axis.title.x = element_text(family = "Arial",
                              size = rel(4),colour = "black"), 
    axis.title.y = element_blank(),
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    legend.text=element_text(size=30),
    legend.title=element_blank()
  )+
  guides(colour = guide_legend(override.aes = list(size=10)))

dev.off()

# #r2
# nullmod <- glm(Success ~ 0, family="binomial", data=mitoVGP)
# 1-logLik(mod.mitoVGP)/logLik(nullmod)
# 
# #dist
# hist(mitoVGP$`Available Pacbio mtDNA reads`)
# hist(sqrt(mitoVGP$`Available Pacbio mtDNA reads`))

# Chi-square test for associtation of assembly success with 
# availability of long reads
Tab2 <- with(mitoVGP,table(Success, I(Available.Pacbio.mtDNA.reads>0)))
chisq.test(Tab2)
# possible alternative tests
fisher.test(Tab2)
chisq.test(Tab2, simulate.p.value = TRUE, B = 10000)

# Chi-square test for the association between tissue type and taxonomic group
chisq.test(out1)
# possible alternative tests
chisq.test(out1, simulate.p.value = TRUE, B = 10000)
# fisher.test(out1) does not work

# keep only categories with N > 3
filtered_mitoVGP <- mitoVGP %>% group_by(Tissue.type.Pacbio) %>% filter(n() > 3) %>%
                                group_by(Taxonomic.group) %>% filter(n() > 3) %>%                               
                                group_by(DNA.extraction.Pacbio) %>% filter(n() > 3) %>%
                                group_by(Library.prep.fragmentation.Pacbio) %>% filter(n() > 3) %>%
                                group_by(Library.prep.Pacbio) %>% filter(n() > 3)

# turn into factor and remove unused levels
filtered_mitoVGP <- filtered_mitoVGP %>% mutate_if(is.character, as.factor)
filtered_mitoVGP <- filtered_mitoVGP %>% mutate_if(is.factor, droplevels)
# for some reason some columns sometimes are not converted into factor
filtered_mitoVGP$Taxonomic.group <- as.factor(filtered_mitoVGP$Taxonomic.group)
filtered_mitoVGP$Tissue.type.Pacbio <- as.factor(filtered_mitoVGP$Tissue.type.Pacbio)
filtered_mitoVGP$Library.prep.fragmentation.Pacbio <- as.factor(filtered_mitoVGP$Library.prep.fragmentation.Pacbio)

# Chi-square test for the association between tissue type and taxonomic group in the filtered dataset
out1_filt <- with(filtered_mitoVGP, table(Tissue.type.Pacbio, Taxonomic.group))
chisq.test(out1_filt)
# possible alternative tests
chisq.test(out1_filt, simulate.p.value = TRUE, B = 10000)
# fisher.test(out1_filt) does not work

mod.mitoVGP <- lm(Available.Pacbio.mtDNA.reads ~ Tissue.type.Pacbio + Taxonomic.group + 
                    Size.selection.kbp.Pacbio + DNA.extraction.Pacbio + 
                    Library.prep.fragmentation.Pacbio + 
                    Library.prep.Pacbio + Total.raw.data.Gbp, data=filtered_mitoVGP)

# check model assumptions
par(mfrow = c(2, 2))
plot(mod.mitoVGP)

par(mfrow = c(1, 1))
acf(mod.mitoVGP$residuals)

mean(mod.mitoVGP$residuals)
lawstat::runs.test(mod.mitoVGP$residuals)

# summary with amount of explained variance
summary(mod.mitoVGP)
af <- anova(mod.mitoVGP)
af
afss <- af$"Sum Sq"
ml<-cbind(af, PctExp=afss/sum(afss)*100)
ml

# check significance with permutation approach
mod.mitoVGP.perm <- lm(Available.Pacbio.mtDNA.reads ~ Tissue.type.Pacbio + Taxonomic.group + 
                         Size.selection.kbp.Pacbio + DNA.extraction.Pacbio + 
                         Library.prep.fragmentation.Pacbio + 
                         Library.prep.Pacbio + Total.raw.data.Gbp, data=filtered_mitoVGP)
summary(mod.mitoVGP.perm)
anova(mod.mitoVGP.perm)


# post-hoc tests
Post.hoc.taxon <- glht(mod.mitoVGP,linfct=mcp(Taxonomic.group = "Tukey"))
summary(Post.hoc.taxon)

Post.hoc.tissue <- glht(mod.mitoVGP, linfct=mcp(Tissue.type.Pacbio = "Tukey"))
summary(Post.hoc.tissue)

Post.hoc.fragmentation <- glht(mod.mitoVGP,linfct=mcp(Library.prep.fragmentation.Pacbio = "Tukey"))
summary(Post.hoc.fragmentation)

# write results
# write.table(ml, file='../Supplementary Tables/ST4.tsv', quote=FALSE, sep='\t')

# 
# ggplot(mod.mitoVGP, aes(Size.selection.kbp.Pacbio, Available.Pacbio.mtDNA.reads)) +
#   geom_point() +
#   stat_smooth(method = lm, se = FALSE, formula = y ~ x) +
#   geom_segment(aes(xend = Size.selection.kbp.Pacbio, yend = .fitted), color = "red", size = 0.3)

# Chi-square test for the association between size selection and success
chisq.test(out4)
# possible alternative tests
chisq.test(out4, simulate.p.value = TRUE, B = 10000)
fisher.test(out4) 
