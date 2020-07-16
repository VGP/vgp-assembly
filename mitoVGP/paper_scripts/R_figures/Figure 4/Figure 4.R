setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(tibble)
library(ggpubr)
library(reshape)
library(taxonomizr)
library(data.table)
library(purrr)
library(forcats)

#fetch NCBI db
prepareDatabase('accessionTaxa.sql')

#load RefSeq data
table1 <- read.csv("mitochondrion AND srcdb_refseq[PROP] AND complete [TITLE] NOT complete cds[TITLE] NOT isolate NOT voucher_03302020.txt", header=FALSE, sep="\t")
names(table1) <- c("Header", "Length", "Accession")
table1$Length <- gsub(" bp circular DNA", "\\1", table1$Length)
table1$Length <- gsub(" bp linear DNA", "\\1", table1$Length)
table1$Length <- as.numeric(gsub(",", "\\1", table1$Length))

table1 = transform(table1, Accession = colsplit(Accession, split = " ", names = c('Accession', 'GI')))
taxaId_Genbank<-accessionToTaxa(as.character(table1$Accession$Accession),"accessionTaxa.sql")
taxons_Genbank<-getTaxonomy(taxaId_Genbank,'accessionTaxa.sql')
table1 <- cbind(table1,taxaId_Genbank, taxons_Genbank)

#load VGP data
table2 <- read.csv("../Figure 4/VGP_lengths.txt", header=TRUE, sep="\t")
taxaId_VGP<-getId(as.character(table2$Species),'accessionTaxa.sql')
taxons_VGP<-getTaxonomy(taxaId_VGP,'accessionTaxa.sql')
table2 <- cbind(table2,taxaId_VGP, taxons_VGP)
table2[table2$ID == "fScaArg1",8]<-"Perciformes"
table2[table2$ID == "fParRan2",8]<-"Perciformes"
table2<-filter(table2, !order %in% c("Forcipulatida","Elopiformes", "Cariamiformes")) %>% droplevels()

#sample RefSeq according to VGP orders
sample_scheme <- data.frame(order = levels(table2$order),
                            n = group_size(table2 %>% group_by(order)))

d1<-bind_rows(replicate(1000, table1 %>% nest(-order) %>% left_join(sample_scheme, by = "order") %>% replace_na(list(n = 0)) %>% mutate(Sample = map2(data, n, sample_n)) %>% unnest(Sample) %>% select(Header, Length), simplify = FALSE))

colnames(d1) <- c("Species","Length")

d2<-data.frame(table2$Species, table2$Length)
colnames(d2) <- c("Species","Length")

#combine datasets
data <- bind_rows("Genbank" = d1, "VGP dataset" = d2, .id = "groups")

#statistics
var.test(d1$Length, d2$Length)

t.result<-t.test(d1$Length, d2$Length, var.equal = TRUE)
t.pvalue<-format.pval(pv = t.result$p.value, digits = 2, eps = 0, nsmall = 3)

ks.result<-ks.test(d1$Length, d2$Length)
ks.pvalue<-format.pval(pv = ks.result$p.value, digits = 2, eps = 0, nsmall = 3)

png("Fig. 4.png", width = 2000, height = 600)

#plot
plot1<-ggplot(data, aes(x = Length)) +
  stat_density(data=d1, aes(x = Length, y = ..density.. * 150 * 200, color="RefSeq dataset",fill="RefSeq dataset"), bw = 150,position="identity",geom="area", size = 2, alpha=0.7) +
  stat_density(data=d2, aes(x = Length, y = ..density.. * 150 * 200, color="VGP dataset", fill="VGP dataset"), bw = 150,position="identity",geom="area", size = 2, alpha=0.5) +
  geom_histogram(data=d2,aes(y = ..count..), 
                 binwidth = 200,  
                 colour = "#dbc446", 
                 fill = "#fae678", alpha=0.5, size=0.6) +
  geom_vline(data=d1, aes(xintercept=mean(Length)),
             color="#335eaa", linetype="dashed", size=2)+
  geom_vline(data=d2, aes(xintercept=mean(Length)),
           color="#ffd700", linetype="dashed", size=2)+
  geom_rug(data=d1, color="#335eaa", size=1.5)+
  geom_rug(data=d2, color="#ffd700", size=1.5)+
  theme(axis.line = element_line(size = 0.6, colour = "black"),
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_line(colour = "black", size=(0.6)),
        axis.ticks.length = unit(.50, "cm"),
        axis.title = element_text(family = "Arial",
                                  size = rel(4)), 
        axis.text = element_text(family = "Arial",
                                 size = rel(4)),
        legend.position = c(0.95, 0.80),
        legend.justification = c("right", "top"),
        legend.text = element_text(family = "Arial",
                                   size = rel(4)),
        legend.key.size = unit(6, 'lines')
        ) +
  labs(y="Counts", x = "Assembly length (bp)") +
  scale_x_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  annotate("text", x = 20000, y = 35, label = paste0("Kolmogorov–Smirnov test: p-value = ",ks.pvalue), family = "Arial",size = 15) +
  scale_color_manual(values = c('VGP dataset' = '#ffd700', 'RefSeq dataset' = '#335eaa')) +
  scale_fill_manual(name="", values = c('VGP dataset' = '#ffd700', 'RefSeq dataset' = '#335eaa')) +
  labs(color="")
plot1
dev.off()
           
