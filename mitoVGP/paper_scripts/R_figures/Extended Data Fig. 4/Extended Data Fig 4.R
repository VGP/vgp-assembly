library(ggplot2)
library(dplyr)
library(plyr)
library(waffle)
library(ggthemes)
library(RColorBrewer)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

mitoVGP <- read.csv("NOVOPlasty_results.tsv", header = TRUE, sep = "\t", na.strings="")

mitoVGP$status <- revalue(mitoVGP$status, c("Out of memory!"="Memory issue"))
mitoVGP$status <- revalue(mitoVGP$status, c("Killed"="Memory issue"))
mitoVGP$status <- revalue(mitoVGP$status, c("Cannot allocate memory"="Memory issue"))

mitoVGP$status <- revalue(mitoVGP$status, c("circular"="Circular assembly"))
mitoVGP$status <- revalue(mitoVGP$status, c("contigs"="Multiple contigs"))

mitoVGP_freq <- mitoVGP %>%
  group_by(status) %>%
  tally()

# Barplot
mitoVGP_plot<- ggplot(mitoVGP_freq, aes(x="", y=n, fill=status))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) +
  xlab("")+ylab("")
mitoVGP_plot

vals <- mitoVGP_freq$n
val_names <- sprintf("%s (%s)", mitoVGP_freq$status, scales::percent(round(vals/sum(vals), 2)))
names(vals) <- val_names

waffle::waffle(vals, colors = c("#DC0000FF", "#00A087FF", "#d5ba4d","white"))

