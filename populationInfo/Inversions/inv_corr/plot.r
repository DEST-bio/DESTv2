DEST=read.table("~/Desktop/inv/DEST.txt",header=T)
DrosEU=read.table("~/Desktop/inv/DrosEU.txt",header=T)
header=colnames(DEST)[2:ncol(DEST)]
header=header[header!="In.3R.K"]
pdf("~/Desktop/inv/inv_corr.pdf",width=12,height=8)
par(mfrow=c(2,3))
for (i in header){
  print(i)
  reg=lm(DEST[[i]]~DrosEU[[i]])
  plot(DrosEU[[i]],DEST[[i]],xlab="DEST",ylab="DrosEU",pch=16,col=rgb(0,0,0,0.2),main=i)
  abline(reg,col="red",lty=2,lwd=2)
  legend("topleft",legend=substitute(expression(italic(R)^2 == myvalue),myvalue =round(summary(reg)$'r.squared',2)),box.col="NA")
}
dev.off()