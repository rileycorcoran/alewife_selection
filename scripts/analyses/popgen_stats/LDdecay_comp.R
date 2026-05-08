args <- commandArgs(trailingOnly = TRUE)

pdf(paste0(out_dir, "/LDdecay_comp.pdf"))
read.table(pdf(paste0(out_dir, "/LDdecay_comp.amos")) -> Eamos;
plot(Eamos[,1]/1000,Eamos[,2], type="l", col="red", main="LD decay", xlab="Distance(Kb)", xlim=c(0,300), ylim=c(0,0.764359875061898), ylab=expression(r^{2}), bty="n", lwd=2)
read.table(pdf(paste0(out_dir, "/LDdecay_comp.bride")) -> Ebride;
lines(Ebride[,1]/1000,Ebride[,2], col="black", lwd=2)
read.table(pdf(paste0(out_dir, "/LDdecay_comp.long")) -> Elong;
lines(Elong[,1]/1000,Elong[,2], col="blue", lwd=2)
read.table(pdf(paste0(out_dir, "/LDdecay_comp.pat")) -> Epat;
lines(Epat[,1]/1000,Epat[,2], col="Purple", lwd=2)
read.table(pdf(paste0(out_dir, "/LDdecay_comp.quon")) -> Equon;
lines(Equon[,1]/1000,Equon[,2], col="green", lwd=2)

legend("topright", c("amos","bride","long","pat","quon"), col=c("red","black","blue","Purple","green"), cex=1, lty=c(1,1,1,1,1), bty="n", lwd=2);
dev.off()


png(pdf(paste0(out_dir, "/LDdecay_comp.png"))
read.table(pdf(paste0(out_dir, "/LDdecay_comp.amos")) -> Eamos;
plot(Eamos[,1]/1000,Eamos[,2], type="l", col="red", main="LD decay", xlab="Distance(Kb)", xlim=c(0,300), ylim=c(0,0.764359875061898), ylab=expression(r^{2}), bty="n", lwd=2)
read.table(pdf(paste0(out_dir, "/LDdecay_comp.bride")) -> Ebride;
lines(Ebride[,1]/1000,Ebride[,2], col="black", lwd=2)
read.table(pdf(paste0(out_dir, "/LDdecay_comp.long")) -> Elong;
lines(Elong[,1]/1000,Elong[,2], col="blue", lwd=2)
read.table(pdf(paste0(out_dir, "/LDdecay_comp.pat")) -> Epat;
lines(Epat[,1]/1000,Epat[,2], col="Purple", lwd=2)
read.table(pdf(paste0(out_dir, "/LDdecay_comp.quon")) -> Equon;
lines(Equon[,1]/1000,Equon[,2], col="green", lwd=2)

legend("topright", c("amos","bride","long","pat","quon"), col=c("red","black","blue","Purple","green"), cex=1, lty=c(1,1,1,1,1), bty="n", lwd=2);
dev.off()

