## SEQ CLEANING
Phyx        <- log10(c(0.01, 0.03, 0.31, 3.42, 37.12));
Gblocks     <- log10(c(0.05, 0.30, 3.22, 53.77, 814.39));
Phyutility  <- log10(c(0.65, 0.76, 8.99,114.84));
x           <- 1:5;
phyutilityx <- log10(c(10,100,1000,10000));
plot(x, Phyx, xlab="log10(sequences)",ylab="log10(time)",type="o", col="blue",
     ylim=c(-2,5), main="Sequence Clean Comparison");

lines(x,Gblocks, type="o", col="orange");
lines(phyutilityx,Phyutility, type="o", col="green");
legend("topleft", c("pxclnsq","Gblocks", "phyutility"), cex=1, 
       col=c("blue", "orange", "green"), pch=21, lty=1, bty="n");

## AA -> CDN
Phyx    <- log10(c(0.01, 0.01, 0.05, 0.52, 5.66));
PAL2NAL <- log10(c(0.05, 0.30, 3.06, 33.39, 369.54));
x       <- 1:5;
plot(x,Phyx, xlab="log10(sequences)",ylab="log10(time)",type="o", col="blue",
     ylim=c(-2,3), main="AA to CDN");
lines(x,PAL2NAL, type="o", col="orange");
legend("topleft", c("pxaatocdn", "PAL2NAL"), cex=1, 
       col=c("blue","orange"), pch=21, lty=1, bty="n");

## MCMC LOGS
x     <- 1:5;
phyx  <- c(4, 8, 11, 13, 16);
lcom  <- c((9*60)+20, (16*60)+45, (25*60)+11, (32*60)+44, (38*60)+47);
lcom2 <- c((13*60)+8, (27*60)+3, (41*60)+45, (54*60)+20, (65*60)+47);
plot(x, phyx, xlab="Number of Log Files",ylab="Time (seconds)",type="o", col="blue",
     ylim=c(0,3950), main="MCMC Log Manipulation Comparison");
lines(x, lcom, type="o", col="orange")
lines(x, lcom2, type="o", col="green")
legend("topleft", c("pxlog", "logcombiner v1.8.2", "logcombiner v2.3.2"), cex=1, 
       col=c("blue", "orange", "green"), pch=21, lty=1, bty="n");

