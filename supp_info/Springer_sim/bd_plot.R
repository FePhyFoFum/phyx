
############################################################
## Figure 1 (bd sim)
############################################################

aa <- as.numeric(read.table("springer_sim_ages.txt")[,1]);
nn <- as.numeric(read.table("springer_sim_sizes.txt")[,1]);

par(mfrow=c(1,2));
hist(aa, xlim=rev(c(0,range(aa)[2])), axes=FALSE, ylab="", xlab="Age of root",
     prob=TRUE, main="", ylim=c(0,0.04));
lines(density(aa, adjust=2), lty="dashed", col="red", lwd=1);
axis(side=1);
axis(side=4);

hist(nn, xlab="Number of extant terminals", prob=TRUE,  main="",
     ylim=c(0,0.0013));
lines(density(nn, adjust=2), lty="dashed", col="red", lwd=1);

## reduce margins ##
par(mar=c(4,1,1,3)+0.1)
hist(aa, xlim=rev(c(0,range(aa)[2])), axes=FALSE, ylab="", xlab="Age of root",
     prob=TRUE, main="", ylim=c(0,0.04));
lines(density(aa, adjust=2), lty="dashed", col="red", lwd=1);
axis(side=1);
axis(side=4);

par(mar=c(4,3,1,1)+0.1)
hist(nn, xlab="Number of extant terminals", ylab="", prob=TRUE, main="",
     ylim=c(0,0.0013));
lines(density(nn, adjust=2), lty="dashed", col="red", lwd=1);

## For comparison with empirical:
require(ape);
phy <- read.tree("springer.tre");
