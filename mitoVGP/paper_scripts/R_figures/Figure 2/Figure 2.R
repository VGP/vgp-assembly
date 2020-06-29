library(shape)
library(randomcoloR)
library(DescTools)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

linMap <- function(x)
  (x*3/max(x))

colorize <- function(x)
{
  if(x == "A"){return("#007FFF")}
  else if(x == "T"){return("#FFF700")}
  else if(x == "G"){return("green")}
  else if(x == "C"){return("#960018")}
}

############rLacAgi1############

png("Fig. 2a.png", width = 4500, height = 5000) 

xlim <- c(0 , 6100)
ylim <- c(0, 60)
plot(0, type = "n", xlim = xlim, ylim = ylim,
     main = "", axes=FALSE, xlab = "", ylab = "")

lwd <- 150

arr.width <- 4.5
arr.length <- 1.4

##loading data

rLacAgi1_pacbio<-read.csv("rLacAgi1_pacbio.txt", na.strings = c(".", "N"), sep="\t", header = FALSE)
names(rLacAgi1_pacbio) <- c("dataset", "start","end","type","name","strand")

vcolor_easy = c("#d5ba4d", "#4DBBD5FF", "#DC0000FF", "#00A087FF")[rLacAgi1_pacbio$type]
varr.length = c(arr.length, arr.length, arr.length, arr.length)[rLacAgi1_pacbio$type]
varr.type = c("triangle", "triangle", "triangle","triangle")[rLacAgi1_pacbio$type]

vx0 <- rLacAgi1_pacbio$start
vx1 <- rLacAgi1_pacbio$end

rLacAgi1_refseq<-read.csv("rLacAgi1_refseq.txt", na.strings = c(".", "N"), sep="\t", header = FALSE)
names(rLacAgi1_refseq) <- c("dataset", "start","end","type","name","strand")

rcolor_easy = c("black","#d5ba4d", "#4DBBD5FF", "#DC0000FF", "#00A087FF")[rLacAgi1_refseq$type]
rarr.length = c(0, arr.length, arr.length, arr.length, arr.length)[rLacAgi1_refseq$type]
rarr.type = c("NA","triangle","triangle", "triangle", "triangle")[rLacAgi1_refseq$type]

rx0 <- rLacAgi1_refseq$start
rx1 <- rLacAgi1_refseq$end

###RECT###
#1
vy0 <- rep(30,13)
vy1 <- rep(30,13)
ry0 <- rep(40,13)
ry1 <- rep(40,13)

rect(vx0[1], vy0[1], rx1[8], ry1[8], density = NULL, angle = 45,
     col = rgb(255, 165, 0, max = 255, alpha = 50), border = NA, lty = par("lty"), lwd = par("lwd"))
#2
vy0 <- rep(20,13)
vy1 <- rep(20,13)
ry0 <- rep(30,13)
ry1 <- rep(30,13)

rect(vx0[1], vy0[1], rx1[5], ry1[5], density = NULL, angle = 45,
     col = rgb(255, 165, 0, max = 255, alpha = 50), border = NA, lty = par("lty"), lwd = par("lwd"))

polygon(c(rx1[5],rx1[2]-207,vx1[5],vx1[9]), c(vy1[1],vy1[1],ry1[1],ry1[1]), density = NULL, angle = 45,
        col = rgb(255, 165, 0, max = 255, alpha = 50), border = NA)

rect(vx0[10], vy0[10], rx1[8], ry1[8], density = NULL, angle = 45,
     col = rgb(255, 165, 0, max = 255, alpha = 50), border = NA, lty = par("lty"), lwd = par("lwd"))

#3
vy0 <- rep(10,13)
vy1 <- rep(10,13)
ry0 <- rep(20,13)
ry1 <- rep(20,13)

rect(vx0[1], vy0[1], rx1[5], ry1[5], density = NULL, angle = 45,
     col = rgb(255, 165, 0, max = 255, alpha = 50), border = NA, lty = par("lty"), lwd = par("lwd"))

rect(vx0[10], vy0[10], rx1[8], ry1[8], density = NULL, angle = 45,
     col = rgb(255, 165, 0, max = 255, alpha = 50), border = NA, lty = par("lty"), lwd = par("lwd"))

####REF####

ry0 <- rep(10,13)
ry1 <- rep(10,13)

for(i in 1:length(rx0)) {
  if(rLacAgi1_refseq$strand[i] == "+"){
    if(!(rLacAgi1_refseq$type[i] == "tRNA")){
      segments(rx0[i], ry0[i], rx1[i]-70, ry1[i], col = rcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(rx0[i], ry0[i], rx1[i], ry1[i], arr.length = rarr.length[i], arr.width = arr.width, code = 2,
           arr.col = rcolor_easy[i], lcol = rcolor_easy[i], segment = FALSE, arr.type = rarr.type[i], arr.adj = 1)
    
  }else if (rLacAgi1_refseq$strand[i] == "-"){
    if(!(rLacAgi1_refseq$type[i] == "tRNA")){
      segments(rx1[i], ry0[i], rx0[i]+70, ry1[i], col = rcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(rx1[i], ry0[i], rx0[i], ry1[i], arr.length = rarr.length[i], arr.width = arr.width, code = 2,
           arr.col = rcolor_easy[i], lcol = rcolor_easy[i], segment = FALSE, arr.type = rarr.type[i], arr.adj = 1)
  }else{
    
    segments(rx1[i], ry0[i], rx0[i]+70, ry1[i], col = rcolor_easy[i], lend=1, lwd=10, lty=2)
    
  }
}

#y<-read.csv("rLacAgi1_RefSeq.snp", sep="\t", header = FALSE)

#for(i in 1:length(y[,1])) {
  
#  segments(y[i,1], 9.25+y[i,4], y[i,1], 10.75-y[i,4], col = colorize(y[i,3]), lwd=1)
  
#}

repeat1<-cbind(seq(from = 1376, to = 1376+199, by = 1))
n=100

color = ""

for(i in 1:nrow(repeat1)) {
  color[i]<-randomColor(hue="blue")
}

DrawArc(x = repeat1+50, y = rep(11-2,length(repeat1)), rx = rep(40,length(repeat1)), ry = -rep(4,length(repeat1)),
        theta.1 = 0.145, theta.2 = 3, nv = 1000,
        col = color, lty = par("lty"), lwd = par("lwd"),
        plot = TRUE)

####NOVOPlasty####

ry0 <- rep(20,13)
ry1 <- rep(20,13)

for(i in 1:length(rx0)) {
  if(rLacAgi1_refseq$strand[i] == "+"){
    if(!(rLacAgi1_refseq$type[i] == "tRNA")){
      segments(rx0[i], ry0[i], rx1[i]-70, ry1[i], col = rcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(rx0[i], ry0[i], rx1[i], ry1[i], arr.length = rarr.length[i], arr.width = arr.width, code = 2,
           arr.col = rcolor_easy[i], lcol = rcolor_easy[i], segment = FALSE, arr.type = rarr.type[i], arr.adj = 1)
    
  }else if (rLacAgi1_refseq$strand[i] == "-"){
    if(!(rLacAgi1_refseq$type[i] == "tRNA")){
      segments(rx1[i], ry0[i], rx0[i]+70, ry1[i], col = rcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(rx1[i], ry0[i], rx0[i], ry1[i], arr.length = rarr.length[i], arr.width = arr.width, code = 2,
           arr.col = rcolor_easy[i], lcol = rcolor_easy[i], segment = FALSE, arr.type = rarr.type[i], arr.adj = 1)
  }else{
    
    segments(rx1[i], ry0[i], rx0[i]+70, ry1[i], col = rcolor_easy[i], lend=1, lwd=10, lty=2)
    
  }
}

repeat1<-cbind(seq(from = 1376, to = 1376+199, by = 1))
n=100

color = ""

for(i in 1:nrow(repeat1)) {
  color[i]<-randomColor(hue="blue")
}

DrawArc(x = repeat1+50, y = rep(11+10-2,length(repeat1)), rx = rep(40,length(repeat1)), ry = -rep(4,length(repeat1)),
        theta.1 = 0.145, theta.2 = 3, nv = 1000,
        col = color, lty = par("lty"), lwd = par("lwd"),
        plot = TRUE)

####VGP_nanopore####

vy0 <- rep(30,13)
vy1 <- rep(30,13)

for(i in 1:length(vx0)) {
  if(rLacAgi1_pacbio$strand[i] == "+"){
    if(!(rLacAgi1_pacbio$type[i] == "tRNA")){
      segments(vx0[i], vy0[i], vx1[i]-70, vy1[i], col = vcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(vx0[i], vy0[i], vx1[i], vy1[i], arr.length = varr.length[i], arr.width = arr.width, code = 2,
           arr.col = vcolor_easy[i], lcol = vcolor_easy[i], segment = FALSE, arr.type = varr.type[i], arr.adj = 1)
    
  }else{
    if(!(rLacAgi1_pacbio$type[i] == "tRNA")){
      segments(vx1[i], vy0[i], vx0[i]+70, vy1[i], col = vcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(vx1[i], vy0[i], vx0[i], vy1[i], arr.length = varr.length[i], arr.width = arr.width, code = 2,
           arr.col = vcolor_easy[i], lcol = vcolor_easy[i], segment = FALSE, arr.type = varr.type[i], arr.adj = 1)
  }
  
}

y<-read.csv("rLacAgi1_ont.cov", sep="\t", header = FALSE)
y<-linMap(y)

for(i in 1:tail(rx1, n=1)) {
  
  segments(i, 31.5, i+1, 31.5+y[i,], col = "gray")
  
}

repeat1<-cbind(seq(from = 1376, to = 1376+199, by = 1))
n=100

color = ""

for(i in 1:nrow(repeat1)) {
  color[i]<-randomColor(hue="blue")
}

DrawArc(x = repeat1+50, y = rep(11+20-2,length(repeat1)), rx = rep(40,length(repeat1)), ry = -rep(4,length(repeat1)),
        theta.1 = 0.145, theta.2 = 3, nv = 1000,
        col = color, lty = par("lty"), lwd = par("lwd"),
        plot = TRUE)

repeat2<-cbind(seq(from = 3379, to = 3577, by = 1))
n=100

color = ""

for(i in 1:nrow(repeat1)) {
  color[i]<-randomColor(hue="blue")
}

DrawArc(x = repeat2+50, y = rep(11+20-2,length(repeat1)), rx = rep(40,length(repeat1)), ry = -rep(4,length(repeat1)),
        theta.1 = 0.145, theta.2 = 3, nv = 1000,
        col = color, lty = par("lty"), lwd = par("lwd"),
        plot = TRUE)

####VGP_pacbio####

vy0 <- rep(40,13)
vy1 <- rep(40,13)

for(i in 1:length(vx0)) {
  if(rLacAgi1_pacbio$strand[i] == "+"){
    if(!(rLacAgi1_pacbio$type[i] == "tRNA")){
      segments(vx0[i], vy0[i], vx1[i]-70, vy1[i], col = vcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(vx0[i], vy0[i], vx1[i], vy1[i], arr.length = varr.length[i], arr.width = arr.width, code = 2,
           arr.col = vcolor_easy[i], lcol = vcolor_easy[i], segment = FALSE, arr.type = varr.type[i], arr.adj = 1)
    
  }else{
    if(!(rLacAgi1_pacbio$type[i] == "tRNA")){
      segments(vx1[i], vy0[i], vx0[i]+70, vy1[i], col = vcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(vx1[i], vy0[i], vx0[i], vy1[i], arr.length = varr.length[i], arr.width = arr.width, code = 2,
           arr.col = vcolor_easy[i], lcol = vcolor_easy[i], segment = FALSE, arr.type = varr.type[i], arr.adj = 1)
  }
  
}

y<-read.csv("rLacAgi1_pacbio.cov", sep="\t", header = FALSE)
y<-linMap(y)

for(i in 1:tail(rx1, n=1)) {
  
  segments(i, 41.5, i+1, 41.5+y[i,], col = "gray")
  
}

repeat1<-cbind(seq(from = 1376, to = 1376+199, by = 1))
n=100

color = ""

for(i in 1:nrow(repeat1)) {
  color[i]<-randomColor(hue="blue")
}

DrawArc(x = repeat1+50, y = rep(11+30-2,length(repeat1)), rx = rep(40,length(repeat1)), ry = -rep(4,length(repeat1)),
        theta.1 = 0.145, theta.2 = 3, nv = 1000,
        col = color, lty = par("lty"), lwd = par("lwd"),
        plot = TRUE)

repeat2<-cbind(seq(from = 3379, to = 3577, by = 1))
n=100

color = ""

for(i in 1:nrow(repeat1)) {
  color[i]<-randomColor(hue="blue")
}

DrawArc(x = repeat2+50, y = rep(11+30-2,length(repeat1)), rx = rep(40,length(repeat1)), ry = -rep(4,length(repeat1)),
        theta.1 = 0.145, theta.2 = 3, nv = 1000,
        col = color, lty = par("lty"), lwd = par("lwd"),
        plot = TRUE)

###LABS###

segments(0, 45, 6081, 45, col = "black", lwd=6)

text(labels="mitoVGP", x=150, y=37.5, cex=10, font=1, family="Arial", pos = 4)
text(labels="mitoVGP", x=150, y=27.5, cex=10, font=1, family="Arial", pos = 4)
text(labels="(Pacbio)", x=150, y=35.5, cex=10, font=1, family="Arial", pos = 4)
text(labels="(Nanopore)", x=150, y=25.5, cex=10, font=1, family="Arial", pos = 4)
text(labels="NOVOPlasty", x=150, y=17.5, cex=10, font=1, family="Arial", pos = 4)
text(labels="NC_021766.1", x=150, y=7.5, cex=10, font=1, family="Arial", pos = 4)

text(labels="a", x=20, y=48.7, cex=13, font=2, family="Arial", pos = 4)
text(labels="L. agilis", x=250, y=48.65, cex=13, font=3, family="Arial", pos = 4)

start = 15000
offset = 15000-14032

labels1=seq(from = start/1000, to = 18, by = 1)
x1=seq(from = 0, to = 4062, by = 1000)
x2=seq(from = start-800, to = start + 5062-1000, by = 200)

for(i in 1:length(x1)) {
  
  text(labels=labels1[i], x=x1[i]+offset, y=47, cex=10, font=1, family="Arial", pos = 1)  
  segments(x1[i]+offset, 45, x1[i]+offset, 45.5, col = "black", lwd=6)
  
}

for(i in x2) {
  
  segments(i-start+offset, 45, i-start+offset, 45.25, col = "black", lwd=6)
  
}

labels2=seq(from = 0, to = 2,by = 1)
x3=seq(from = 5063, to = 6081, by = 1000)
x4=seq(from = 5063, to = 6081, by = 200)

for(i in 1:length(x3)) {
  
  text(labels=labels2[i], x=x3[i], y=47, cex=10, font=1, family="Arial", pos = 1)  
  segments(x3[i], 45, x3[i], 45.5, col = "black", lwd=6)
  
}

for(i in x4) {
  
  segments(i, 45, i, 45.25, col = "black", lwd=6)
  
}

dev.off()

###########################



###########################

png("Fig. 2b.png", width = 4500, height = 2500) 

###########################

xlim <- c(0 , 6500)
ylim <- c(0, 30)
plot(0, type = "n", xlim = xlim, ylim = ylim,
     main = "", axes=FALSE, ylab='', xlab='')

text(labels="b", x=20, y=28.8, cex=13, font=2, family="Arial")
text(labels="A. chrysaetos", x=250, y=28.6, cex=13, font=3, family="Arial", pos = 4)

lwd <- 150

arr.width <- 4.5
arr.length <- 1.25

############bAquChr1############

##loading data

bAquChr1<-read.csv("bAquChr1.txt", na.strings = c(".", "N"), sep="\t", header = FALSE)
names(bAquChr1) <- c("dataset", "start","end","type","name","stat","strand")

vcolor_easy = c("#d5ba4d", "#4DBBD5FF", "#DC0000FF", "#00A087FF")[bAquChr1$type]
varr.length = c(arr.length, arr.length, arr.length, arr.length)[bAquChr1$type]
varr.type = c("triangle", "triangle", "triangle", "triangle")[bAquChr1$type]

vx0 <- bAquChr1$start
vy0 <- rep(20,13)
vx1 <- bAquChr1$end
vy1 <- rep(20,13)

NC_024087<-read.csv("NC_024087.1.txt", na.strings = c(".", "N"), sep="\t", header = FALSE)
names(NC_024087) <- c("dataset", "start","end","type","name","stat","strand")

rcolor_easy = c("black", "#d5ba4d", "#4DBBD5FF", "#DC0000FF", "#00A087FF")[NC_024087$type]
rarr.length = c(0, arr.length, arr.length, arr.length, arr.length)[NC_024087$type]
rarr.type = c("NA","triangle", "triangle", "triangle", "triangle")[NC_024087$type]

rx0 <- NC_024087$start
ry0 <- rep(10,14)
rx1 <- NC_024087$end
ry1 <- rep(10,14)

###LABS###

segments(0, 25, 6339, 25, col = "black", lwd=6)

text(labels="NC_024087.1", x=150, y=7.5, cex=10, font=1, family="Arial", pos = 4)
text(labels="Pacbio", x=150, y=17.5, cex=10, font=1, family="Arial", pos = 4)

start = 15000
offset = 15000-14529

labels1=seq(from = start/1000, to = 18, by = 1)
x1=seq(from = 0, to = 3733, by = 1000)
x2=seq(from = start-400, to = start + 3733-400, by = 200)

for(i in 1:length(x1)) {
  
  text(labels=labels1[i], x=x1[i]+offset, y=27, cex=10, font=1, family="Arial", pos = 1)  
  segments(x1[i]+offset, 25, x1[i]+offset, 25.5, col = "black", lwd=6)
  
}

for(i in x2) {
  
  segments(i-start+offset, 25, i-start+offset, 25.25, col = "black", lwd=6)
  
}

labels2=seq(from = 0, to = 2,by = 1)
x3=seq(from = 3734, to = 6339, by = 1000)
x4=seq(from = 3734, to = 6339, by = 200)

for(i in 1:length(x3)) {
  
  text(labels=labels2[i], x=x3[i], y=27, cex=10, font=1, family="Arial", pos = 1)  
  segments(x3[i], 25, x3[i], 25.5, col = "black", lwd=6)
  
}

for(i in x4) {
  
  segments(i, 25, i, 25.25, col = "black", lwd=6)
  
}

###RECT###

rect(vx0[8]+10, vy0[8], rx1[12], ry1[12], density = NULL, angle = 45,
     col = rgb(255, 165, 0, max = 255, alpha = 50), border = NA, lty = par("lty"), lwd = par("lwd"))

rect(vx0[1], vy0[1], rx1[7], ry1[8], density = NULL, angle = 45,
     col = rgb(255, 165, 0, max = 255, alpha = 50), border = NA, lty = par("lty"), lwd = par("lwd"))

####REF####

for(i in 1:length(rx0)) {
  if(NC_024087$strand[i] == "+"){
    if(!(NC_024087$type[i] == "tRNA")){
      segments(rx0[i], ry0[i], rx1[i]-60, ry1[i], col = rcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(rx0[i], ry0[i], rx1[i], ry1[i], arr.length = rarr.length[i], arr.width = arr.width, code = 2,
           arr.col = rcolor_easy[i], lcol = rcolor_easy[i], segment = FALSE, arr.type = rarr.type[i], arr.adj = 1)
    
  }else if (NC_024087$strand[i] == "-"){
    if(!(NC_024087$type[i] == "tRNA")){
      segments(rx1[i], ry0[i], rx0[i]+60, ry1[i], col = rcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(rx1[i], ry0[i], rx0[i], ry1[i], arr.length = rarr.length[i], arr.width = arr.width, code = 2,
           arr.col = rcolor_easy[i], lcol = rcolor_easy[i], segment = FALSE, arr.type = rarr.type[i], arr.adj = 1)
  }else{
    
    segments(rx1[i], ry0[i], rx0[i]+60, ry1[i], col = rcolor_easy[i], lend=1, lwd=10, lty=2)
    
  }
  
  
  
}

repeat1<-cbind(seq(from = 2472-98+70, to = 2472+70, by = 1))
n=100

for(i in 1:nrow(repeat1)) {
  color[i]<-randomColor(hue="blue")
}

DrawArc(x = repeat1-110, y = rep(11-2,length(repeat1)), rx = rep(40,length(repeat1)), ry = -rep(4,length(repeat1)),
        theta.1 = 0.145, theta.2 = 3, nv = 1000,
        col = color, lty = par("lty"), lwd = par("lwd"),
        plot = TRUE)

####VGP####

for(i in 1:length(vx0)) {
  if(bAquChr1$strand[i] == "+"){
    if(!(bAquChr1$type[i] == "tRNA")){
      segments(vx0[i], vy0[i], vx1[i]-60, vy1[i], col = vcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(vx0[i], vy0[i], vx1[i], vy1[i], arr.length = varr.length[i], arr.width = arr.width, code = 2,
           arr.col = vcolor_easy[i], lcol = vcolor_easy[i], segment = FALSE, arr.type = varr.type[i], arr.adj = 1)
    
  }else{
    if(!(bAquChr1$type[i] == "tRNA")){
      segments(vx1[i], vy0[i], vx0[i]+60, vy1[i], col = vcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(vx1[i], vy0[i], vx0[i], vy1[i], arr.length = varr.length[i], arr.width = arr.width, code = 2,
           arr.col = vcolor_easy[i], lcol = vcolor_easy[i], segment = FALSE, arr.type = varr.type[i], arr.adj = 1)
  }
  
}

linMap <- function(x)
  (x*3/max(x))

y<-read.csv("bAquChr1.cov", sep="\t", header = FALSE)
y<-linMap(y)

for(i in 1:tail(rx1, n=1)) {
  
  segments(i, 21.5, i+1, 21.5+y[i,], col = "gray")
  
}

repeat1<-cbind(seq(from = 3733-1058+70, to = 3733+70, by = 1))
n=100

color = ""

for(i in 1:nrow(repeat1)) {
  color[i]<-randomColor(hue="blue")
}

DrawArc(x = repeat1-130, y = rep(11+8,length(repeat1)), rx = rep(40,length(repeat1)), ry = -rep(4,length(repeat1)),
        theta.1 = 0.140, theta.2 = 3, nv = 1000,
        col = color, lty = par("lty"), lwd = par("lwd"),
        plot = TRUE)

dev.off() 
###########################

png("Fig. 2c.png", width = 4500, height = 2500) 

###########################

xlim <- c(0 , 6500)
ylim <- c(0, 30)
plot(0, type = "n", xlim = xlim, ylim = ylim,
     main = "", axes=FALSE, ylab='', xlab='')

text(labels="c", x=20, y=28.8, cex=13, font=2, family="Arial")
text(labels="S. habroptilus", x=250, y=28.6, cex=13, font=3, family="Arial", pos = 4)

lwd <- 150

arr.width <- 4.5
arr.length <- 1.25

############bStrHab1############

##loading data

bStrHab1<-read.csv("bStrHab1.txt", na.strings = c(".", "N"), sep="\t", header = FALSE)
names(bStrHab1) <- c("dataset", "start","end","type","name","stat","strand")

vcolor_easy = c("#d5ba4d", "#4DBBD5FF", "#DC0000FF", "#00A087FF")[bStrHab1$type]
varr.length = c(arr.length, arr.length, arr.length, arr.length)[bStrHab1$type]
varr.type = c("triangle", "triangle", "triangle", "triangle")[bStrHab1$type]

vx0 <- bStrHab1$start
vy0 <- rep(20,13)
vx1 <- bStrHab1$end
vy1 <- rep(20,13)

NC005931<-read.csv("NC_005931.1.txt", na.strings = c(".", "N"), sep="\t", header = FALSE)
names(NC005931) <- c("dataset", "start","end","type","name","stat","strand")

rcolor_easy = c("black", "#d5ba4d", "#4DBBD5FF", "#DC0000FF", "#00A087FF")[NC005931$type]
rarr.length = c(0, arr.length, arr.length, arr.length, arr.length)[NC005931$type]
rarr.type = c("NA","triangle", "triangle", "triangle", "triangle")[NC005931$type]

rx0 <- NC005931$start
ry0 <- rep(10,14)
rx1 <- NC005931$end
ry1 <- rep(10,14)

###LABS###

segments(0, 25, 6339, 25, col = "black", lwd=6)

text(labels="NC_005931.1", x=150, y=7.5, cex=10, font=1, family="Arial", pos = 4)
text(labels="Pacbio", x=150, y=17.5, cex=10, font=1, family="Arial", pos = 4)

start = 14000
offset = 14000-13649

labels1=seq(from = start/1000, to = 17, by = 1)
x1=seq(from = 0, to = 3180, by = 1000)
x2=seq(from = start-200, to = start + 4180-200, by = 200)

for(i in 1:length(x1)) {
  
  text(labels=labels1[i], x=x1[i]+offset, y=27, cex=10, font=1, family="Arial", pos = 1)  
  segments(x1[i]+offset, 25, x1[i]+offset, 25.5, col = "black", lwd=6)
  
}

for(i in x2) {
  
  segments(i-start+offset, 25, i-start+offset, 25.25, col = "black", lwd=6)
  
}

labels2=seq(from = 0, to = 2,by = 1)
x3=seq(from = 4181, to = 6339, by = 1000)
x4=seq(from = 4181, to = 6339, by = 200)

for(i in 1:length(x3)) {
  
  text(labels=labels2[i], x=x3[i], y=27, cex=10, font=1, family="Arial", pos = 1)  
  segments(x3[i], 25, x3[i], 25.5, col = "black", lwd=6)
  
}

for(i in x4) {
  
  segments(i, 25, i, 25.25, col = "black", lwd=6)
  
}

###RECT###

rect(vx0[7], vy0[8], rx1[11], ry1[14], density = NULL, angle = 45,
     col = rgb(255, 165, 0, max = 255, alpha = 50), border = NA, lty = par("lty"), lwd = par("lwd"))

rect(vx0[1], vy0[1], rx1[6], ry1[8], density = NULL, angle = 45,
     col = rgb(255, 165, 0, max = 255, alpha = 50), border = NA, lty = par("lty"), lwd = par("lwd"))

####REF####

for(i in 1:length(rx0)) {
  if(NC005931$strand[i] == "+"){
    if(!(NC005931$type[i] == "tRNA")){
      segments(rx0[i], ry0[i], rx1[i]-60, ry1[i], col = rcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(rx0[i], ry0[i], rx1[i], ry1[i], arr.length = rarr.length[i], arr.width = arr.width, code = 2,
           arr.col = rcolor_easy[i], lcol = rcolor_easy[i], segment = FALSE, arr.type = rarr.type[i], arr.adj = 1)
    
  }else if (NC005931$strand[i] == "-"){
    if(!(NC005931$type[i] == "tRNA")){
      segments(rx1[i], ry0[i], rx0[i]+60, ry1[i], col = rcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(rx1[i], ry0[i], rx0[i], ry1[i], arr.length = rarr.length[i], arr.width = arr.width, code = 2,
           arr.col = rcolor_easy[i], lcol = rcolor_easy[i], segment = FALSE, arr.type = rarr.type[i], arr.adj = 1)
  }else{
    
    segments(rx1[i], ry0[i], rx0[i]+60, ry1[i], col = rcolor_easy[i], lend=1, lwd=10, lty=2)
    
  }
  
  
  
}

####VGP####

for(i in 1:length(vx0)) {
  if(bStrHab1$strand[i] == "+"){
    if(!(bStrHab1$type[i] == "tRNA")){
      segments(vx0[i], vy0[i], vx1[i]-60, vy1[i], col = vcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(vx0[i], vy0[i], vx1[i], vy1[i], arr.length = varr.length[i], arr.width = arr.width, code = 2,
           arr.col = vcolor_easy[i], lcol = vcolor_easy[i], segment = FALSE, arr.type = varr.type[i], arr.adj = 1)
    
  }else{
    if(!(bStrHab1$type[i] == "tRNA")){
      segments(vx1[i], vy0[i], vx0[i]+60, vy1[i], col = vcolor_easy[i], lend=1, lwd=lwd)
    }
    Arrows(vx1[i], vy0[i], vx0[i], vy1[i], arr.length = varr.length[i], arr.width = arr.width, code = 2,
           arr.col = vcolor_easy[i], lcol = vcolor_easy[i], segment = FALSE, arr.type = varr.type[i], arr.adj = 1)
  }
  
}

linMap <- function(x)
  (x*3/max(x))

y<-read.csv("bStrHab1.cov", sep="\t", header = FALSE)
y<-linMap(y)

for(i in 1:tail(rx1, n=1)) {
  
  segments(i, 21.5, i+1, 21.5+y[i,], col = "gray")
  
}

repeat1<-cbind(seq(from = 3074, to = 3074+985, by = 1))
n=100

color = ""

for(i in 1:nrow(repeat1)) {
  color[i]<-randomColor(hue="blue")
}

DrawArc(x = repeat1-130, y = rep(11+8,length(repeat1)), rx = rep(40,length(repeat1)), ry = -rep(4,length(repeat1)),
        theta.1 = 0.140, theta.2 = 3, nv = 1000,
        col = color, lty = par("lty"), lwd = par("lwd"),
        plot = TRUE)

dev.off() 
###########################




###########################

library(circlize)
library("scales")
library(randomcoloR)

setwd("~/Documents/VGP/com/papers/mitoVGP/submission_version/Figure 2/VGP_noref/")

############bCicMag1############

##loading data

bCicMag1<-read.csv("bCicMag1/circos.tsv", na.strings = c(".", "N"), sep="\t", header = FALSE)
names(bCicMag1) <- c("dataset", "type","name","met","start","end")

color_easy = c("#d5ba4d", "#4DBBD5FF", "#DC0000FF", "#00A087FF","#3C5488FF")[bCicMag1$type]

sectors = cbind(bCicMag1$start, bCicMag1$end)

tbl<-data.frame(1:21477)
colnames(tbl) <- "BP"

dat<-do.call(rbind, apply(tbl, 1, function(x) {
  if(length(idx <- which(x["BP"] >= bCicMag1$start & x["BP"] < bCicMag1$end)) > 0) {
    cbind(rbind(x), bCicMag1[idx,], row.names = NULL)
  } else {cbind(rbind(x), bCicMag1[1,][NA,], row.names = NULL)}
}))

cov<-scan("bCicMag1/cov.txt", numeric(), quote = "")

data<-cbind(dat,cov)

circos.clear()

png("Fig. 2d.png", width = 2000, height = 2000) 

##first circle

circos.par(start.degree = 90, cell.padding = c(0.02, 0, 0.02, 0), gap.degree=0, canvas.xlim=c(-1.1, 1.1))
circos.initialize("scale", xlim = c(1,21477))
circos.track(ylim = c(0, 2),bg.border = NA)

labels<-seq(from = 0, to = 21000, by = 1000) 
labels<-lapply(labels, function(x) format(x/1000, big.mark = ",", scientific = FALSE)) 

circos.axis(labels = labels, "top", labels.facing="clockwise", 
            major.at = seq(from = 0, to = 21477, by = 1000), labels.cex = 7, lwd = 10,
            major.tick.length = 0.5)

text(-1, 1, "d", cex=7, font=2, family="Arial")

##second circle

par(new = TRUE) # <- magic
circos.initialize("bCicMag1", xlim = c(1,21477))
circos.track(ylim = c(0, 250),bg.border = NA)
circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))

circos.trackLines(rep("bCicMag1",21477), data$BP, data$cov, area = TRUE, type='s', col="grey", border = NA)

##third circle

par(new = TRUE) # <- magic
circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1.32, 1.32))
circos.initialize(bCicMag1$name, xlim = sectors)

circos.track(ylim = c(1, 2), panel.fun = function(x, y, ...) {
  circos.arrow(CELL_META$xlim[1], CELL_META$xlim[2], 
               arrow.head.width = CELL_META$yrange*0.75, arrow.head.length = ux(0.45, "cm"),
               col = color_easy[CELL_META$sector.numeric.index])
}, bg.border = NA, "track.height" = 0.28)

circos.link("cob_1", c(14266,16788), "cob_2", c(17515,20020), col = rgb(255, 165, 0, max = 255, alpha = 50),
            border = rgb(255, 165, 0, max = 255, alpha = 50), h = 1)

repeat1<-data.frame(cbind(seq(from = 16789, to = 17508, by = 15),seq(from = 16789, to = 17508, by = 15)+1,
                          seq(from = 20021, to = 20544, by = 15),seq(from = 20021, to = 20544, by = 15)+1
))

for(i in 1:nrow(repeat1)) {
  color<-randomColor(hue="blue")
  circos.link("cob_1", cbind(repeat1[i,]$X1,repeat1[i,]$X2), "cob_2", cbind(repeat1[i,]$X3,repeat1[i,]$X4),
              col = color, border = color, h = 1)
}

repeat2<-data.frame(cbind(seq(from = 20635, to = 21477, by = 15),seq(from = 20635, to = 21477, by = 15)+1,
                          seq(from = 20635, to = 21477, by = 15),seq(from = 20635, to = 21477, by = 15)+1
))

for(i in 1:nrow(repeat2)) {
  color<-randomColor(hue="pink")
  circos.link("cob_1", cbind(repeat2[i,]$X1,repeat2[i,]$X2), "cob_2", cbind(repeat2[i,]$X3,repeat2[i,]$X4),
              col = color, border = color, h = 0.2)
}

text(0, 0, "C. maguari", cex = 10, font=3)

dev.off()

#############################################

############fAntMac1############

##loading data

fAntMac1<-read.csv("fAntMac1/circos.tsv", na.strings = c(".", "N"), sep="\t", header = FALSE)
names(fAntMac1) <- c("dataset", "type","name","met","start","end")

color_easy = c("#d5ba4d", "#4DBBD5FF", "#DC0000FF", "#00A087FF","#3C5488FF")[fAntMac1$type]

sectors = cbind(fAntMac1$start, fAntMac1$end)

tbl<-data.frame(1:19188)
colnames(tbl) <- "BP"

dat<-do.call(rbind, apply(tbl, 1, function(x) {
  if(length(idx <- which(x["BP"] >= fAntMac1$start & x["BP"] < fAntMac1$end)) > 0) {
    cbind(rbind(x), fAntMac1[idx,], row.names = NULL)
  } else {cbind(rbind(x), fAntMac1[1,][NA,], row.names = NULL)}
}))

cov<-scan("fAntMac1/cov.txt", numeric(), quote = "")

data<-cbind(dat,cov)

circos.clear()

png("Fig. 2e.png", width = 2000, height = 2000) 

##first circle

circos.par(start.degree = 90, cell.padding = c(0.02, 0, 0.02, 0), gap.degree=0, canvas.xlim=c(-1.1, 1.1))
circos.initialize("scale", xlim = c(1,19188))
circos.track(ylim = c(0, 2),bg.border = NA)

labels<-seq(from = 0, to = 19188, by = 1000) 
labels<-lapply(labels, function(x) format(x/1000, big.mark = ",", scientific = FALSE)) 

circos.axis(labels = labels, "top", labels.facing="clockwise", 
            major.at = seq(from = 0, to = 19188, by = 1000), labels.cex = 7, lwd = 10,
            major.tick.length = 0.5)

text(-1, 1, "e", cex=7, font=2, family="Arial")

##second circle

par(new = TRUE) # <- magic
circos.initialize("fAntMac1", xlim = c(1,19188))
circos.track(ylim = c(0, 30),bg.border = NA)
circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))

circos.trackLines(rep("fAntMac1",19188), data$BP, data$cov, area = TRUE, type='s', col="grey", border = NA)

##third circle

par(new = TRUE) # <- magic
circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1.32, 1.32))
circos.initialize(fAntMac1$name, xlim = sectors)

circos.track(ylim = c(1, 1.1), panel.fun = function(x, y, ...) {
  circos.arrow(CELL_META$xlim[1], CELL_META$xlim[2], 
               arrow.head.width = CELL_META$yrange*0.75, arrow.head.length = ux(0.45, "cm"),
               col = color_easy[CELL_META$sector.numeric.index], border = NA)
}, bg.border = NA, "track.height" = 0.28)

circos.link("trnV_1", c(979,1397), "trnV_2", c(1946,2364), col = rgb(255, 165, 0, max = 255, alpha = 50),
            border = rgb(255, 165, 0, max = 255, alpha = 50), h = 20)

repeat1<-data.frame(cbind(seq(from = 16714, to = 17193, by = 15),seq(from = 16714, to = 17193, by = 15)+1,
                          seq(from = 16714, to = 17193, by = 15),seq(from = 16714, to = 17193, by = 15)+1
))

for(i in 1:nrow(repeat1)) {
  color<-randomColor(hue="blue")
  circos.link("OH", cbind(repeat1[i,]$X1,repeat1[i,]$X2), "OH", cbind(repeat1[i,]$X3,repeat1[i,]$X4),
              col = color, border = color, h = 0.2)
}

repeat2<-data.frame(cbind(seq(from = 17787, to = 19135, by = 15),seq(from = 17787, to = 19135, by = 15)+1,
                          seq(from = 17787, to = 19135, by = 15),seq(from = 17787, to = 19135, by = 15)+1
))

for(i in 1:nrow(repeat2)) {
  color<-randomColor(hue="blue")
  circos.link("OH", cbind(repeat2[i,]$X1,repeat2[i,]$X2), "OH", cbind(repeat2[i,]$X3,repeat2[i,]$X4),
              col = color, border = color, h = 0.2)
}

text(0, 0, "A. maculatus", cex=10, font=3)

dev.off() 

################################
