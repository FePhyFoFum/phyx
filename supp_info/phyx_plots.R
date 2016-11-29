## SEQ CLEANING
Phyx        <- log10(c(0.01, 0.03, 0.31, 3.42, 37.12));
Gblocks     <- log10(c(0.05, 0.30, 3.22, 53.77, 814.39));
Phyutility  <- log10(c(0.65, 0.76, 8.99,114.84));
x           <- 1:5;
phyutilityx <- 1:4;
plot(x, Phyx, xlab="log10(sequences)",ylab="log10(seconds)",type="o", col="blue",
     ylim=c(-2, 5), main="");

lines(x,Gblocks, type="o", col="orange");
lines(phyutilityx, Phyutility, type="o", col="green");
legend("topleft", c("pxclsq","Gblocks", "phyutility"), cex=1, 
       col=c("blue", "orange", "green"), pch=21, lty=1, bty="n");

# no log
sPhyx        <- c(0.01, 0.03, 0.31, 3.42, 37.12);
sGblocks     <- c(0.05, 0.30, 3.22, 53.77, 814.39);
sPhyutility  <- c(0.65, 0.76, 8.99,114.84);
x           <- exp(1:5);
phyutilityx <- exp(1:4);
plot(x, sPhyx, xlab="Number Of Sequences)",ylab=" Time(seconds)",type="o", col="blue",
     ylim=c(0, 825), main="");

lines(x, sGblocks, type="o", col="orange");
lines(phyutilityx, sPhyutility, type="o", col="green");
legend("topleft", c("pxclsq","Gblocks", "phyutility"), cex=1, 
       col=c("blue", "orange", "green"), pch=21, lty=1, bty="n");

###############
## AA -> CDN ##
Phyx    <- log10(c(0.01, 0.01, 0.05, 0.52, 5.66));
PAL2NAL <- log10(c(0.05, 0.30, 3.06, 33.39, 369.54));
x       <- 1:5;
plot(x,Phyx, xlab="log10(sequences)",ylab="log10(seconds)",type="o", col="blue",
     ylim=c(-2,3), main="");
lines(x,PAL2NAL, type="o", col="orange");
legend("topleft", c("pxaa2cdn", "PAL2NAL"), cex=1, 
       col=c("blue","orange"), pch=21, lty=1, bty="n");

# no log
sPhyx    <- c(0.01, 0.01, 0.05, 0.52, 5.66);
sPAL2NAL <- c(0.05, 0.30, 3.06, 33.39, 369.54);
x       <- exp(1:5);
plot(x, sPhyx, xlab=" Number Of Sequences)",ylab="Time (seconds)",type="o", col="blue",
     ylim=c(0,375), main="");
lines(x, sPAL2NAL, type="o", col="orange");
legend("topleft", c("pxaa2cdn", "PAL2NAL"), cex=1, 
       col=c("blue","orange"), pch=21, lty=1, bty="n");

###############
## MCMC LOGS ##
x     <- 1:5;
phyx  <- c(4, 8, 11, 13, 16);
lcom  <- c((9*60)+20, (16*60)+45, (25*60)+11, (32*60)+44, (38*60)+47);
lcom2 <- c((13*60)+8, (27*60)+3, (41*60)+45, (54*60)+20, (65*60)+47);
plot(x, phyx, xlab="Number of Log Files",ylab="Time (seconds)",type="o", col="blue",
     ylim=c(0,3950), main="");
lines(x, lcom, type="o", col="orange")
lines(x, lcom2, type="o", col="green")
legend("topleft", c("pxlog", "logcombiner v1.8.2", "logcombiner v2.3.2"), cex=1, 
       col=c("blue", "orange", "green"), pch=21, lty=1, bty="n");

# with log y axis
x     <- 1:5;
lphyx  <- log10(c(4, 8, 11, 13, 16));
llcom  <- log10(c((9*60)+20, (16*60)+45, (25*60)+11, (32*60)+44, (38*60)+47));
llcom2 <- log10(c((13*60)+8, (27*60)+3, (41*60)+45, (54*60)+20, (65*60)+47));
plot(x, lphyx, xlab="Number of Log Files",ylab="log10(seconds)",type="o", col="blue", 
     ylim=c(0,4), main="");
lines(x, llcom, type="o", col="orange")
lines(x, llcom2, type="o", col="green")
legend("topleft", c("pxlog", "logcombiner v1.8.2", "logcombiner v2.3.2"), cex=1, 
       col=c("blue", "orange", "green"), pch=21, lty=1, bty="n");

