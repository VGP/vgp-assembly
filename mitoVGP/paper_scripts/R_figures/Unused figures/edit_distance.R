library(ggplot2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

summary<-read.csv("summary_bCalAnn1.txt", na.strings = c(".", "N"), sep="\t", header = TRUE)

plot(summary$edit_fraction,summary$read_fraction)


p = ggplot() + 
  geom_line(data = summary, aes(x = edit_fraction, y = read_fraction), color = "blue") +
  geom_line(data = summary, aes(x = edit_fraction, y = gbp_fraction), color = "red") +
  xlab('Edit distance/genome length %') +
  ylab('percent.change')

p
