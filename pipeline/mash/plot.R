args<-commandArgs(TRUE)
out_name=paste("mash.", args[[1]], ".png", sep="")

png(out_name, res=300, width=9.6, height=9.6, units="in", pointsize=8, type="cairo", bg="white")
par(mar=c(3.2,3.2,0,0))
par(mgp=c(2.2,1,0))

# get names
key=read.table("key")
labels=key[,1]

# make dendrogram
x = read.table("combined.tbl")
y=x[,2:dim(x)[2]]
z = data.matrix(y)
z[is.infinite(z)]=0
rc = hclust(as.dist(z), method="ward.D2")

# make coloring
x = read.table("combined.tbl");
y=x[,2:dim(x)[2]]
cr3 = data.matrix(y)
cr3=100-(cr3*100)
n=100 # number of steps between 2 colors
mini=min(cr3[])
maxi=max(cr3[row(cr3)!=col(cr3)])
trueMax=max(cr3[])
q25=quantile(cr3[row(cr3)!=col(cr3)],0.25,1)
q50=quantile(cr3[row(cr3)!=col(cr3)],0.5,1)
q75=quantile(cr3[row(cr3)!=col(cr3)],0.75,1)
mini=max(q25-1.5*(q75-q25),0)
maxi=min(q75+1.5*(q75-q25),trueMax)
diff=maxi-mini
palette=colorRampPalette(c("lightyellow", "yellow", "red", "brown", "grey23"))(n = 5*n-1)
breaks=c(seq(mini,mini+diff/4-0.1,length=n), # for lightyellow
                seq(mini+diff/4,mini+diff/2-0.1,length=n), # for yellow
                seq(mini+diff/2,mini+3*diff/4-0.1,length=n), # for red
                seq(mini+3*diff/4,maxi-5,length=n), # for brown
                seq(from=maxi-5+0.1, to=trueMax, length=n))

# plot
library("gplots")
heatmap.2(cr3, Rowv=rev(as.dendrogram(rc)), Colv=as.dendrogram(rc), labRow=as.matrix(labels), labCol=as.matrix(labels), scale="none", distfun=as.dist, col=palette, trace="none", breaks=sort(breaks), margins=c(10,10))
dev.off()
