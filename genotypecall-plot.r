setwd("D:/OneDrive - University College London/A/0ngoing/Heliconius/genotypecall");
GTerror001 = read.table("genotypecallingerror_e0.001.txt", header=T);
GTerror005 = read.table("genotypecallingerror_e0.005.txt", header=T);
GTerror01  = read.table("genotypecallingerror_e0.01.txt", header=T);
GTerror05  = read.table("genotypecallingerror_e0.05.txt", header=T);
GTerror1   = read.table("genotypecallingerror_e0.1.txt", header=T);
par(mfrow=c(1,2))

plot(GTerror001[,1],  GTerror001[,2], type="o", pch=0, col="black", ylim=c(0,0.27), xlab="reads (n)", ylab="genotype calling error")
lines(GTerror005[,1], GTerror005[,2], type="o", pch=1, col="orange")
lines(GTerror01[,1],  GTerror01[,2], type="o", pch=2, col="red")
lines(GTerror05[,1],  GTerror05[,2], type="o", pch=5, col="darkgreen")
lines(GTerror1[,1],   GTerror1[,2], type="o", pch=6, col="blue")
legend("topright", c("e = 0.1", "e = 0.05","e = 0.01","e = 0.005","e = 0.001"), pch=c(6,5,2,1,0), col=c("blue","darkgreen","red","orange","black"))
title("homozygotes")

plot(GTerror001[,1]-0.3,   GTerror001[,3], type="o", pch=0, col="black", ylim=c(0,0.65), xlab="reads (n)", ylab="")
lines(GTerror005[,1]-0.15, GTerror005[,3], type="o", pch=1, col="orange")
lines(GTerror01[,1]+0.00,  GTerror01[,3], type="o", pch=2, col="red")
lines(GTerror05[,1]+0.15,  GTerror05[,3], type="o", pch=5, col="darkgreen")
lines(GTerror1[,1]+0.3,    GTerror1[,3], type="o", pch=6, col="blue")
legend("topright", c("e = 0.1", "e = 0.05","e = 0.01","e = 0.005","e = 0.001"), pch=c(6,5,2,1,0), col=c("blue","darkgreen","red","orange","black"))
title("heterozygotes")
