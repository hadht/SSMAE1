################################################################################
##
## GRAPHS: Adjusting a tSate Space Model to Aedes Aegypti Trap Data    
##
################################################################################
# date: 01/01/2018
#############################################################################

require(chron)
#############################################################################
# set directory:

#############################################################################
### loading:
EstFKVE=read.table("att1EggDengueModelv21.csv",sep=",",dec=".",header=TRUE)
EstFKVE=EstFKVE[,]
EstFKVEpt=read.table("ptt1EggDengueModelv21.csv",sep=",",dec=".",header=TRUE)
EstFKVEpt=EstFKVEpt[,]
EstFKVEytt1=read.table("ytt1EggDengueModelv21.csv",sep=",",dec=".",header=TRUE)
EstFKVEytt1
ytt1=EstFKVEytt1$ytt1
EstFKVEft=read.table("FtEggDengueModelv21.csv",sep=",",dec=".",header=TRUE)
EstFKVEft
Ft=EstFKVEft$Ft
yt=read.table("ytovos-costarica.csv",sep="",dec=",",header=FALSE)
yt=log(yt$V1)
n1=length(yt)

#############################################################################
##
## GRAPHS:
##
############################################################################  
Et=(EstFKVE$Et)
qlils=1.96
liEt=Et-qlils*sqrt(EstFKVEpt$Et)
lsEt=Et+qlils*sqrt(EstFKVEpt$Et)
liEt
Et
lsEt
Etm=matrix(0,n1,3)
Etm[,1]=liEt
Etm[,2]=lsEt
Etm[,3]=Et
Lt=(EstFKVE$lt)
liLt=Lt-qlils*sqrt(EstFKVEpt$lt)
lsLt=Lt+qlils*sqrt(EstFKVEpt$lt)
Ltm=matrix(0,n1,3)
Ltm[,1]=liLt
Ltm[,2]=lsLt
Ltm[,3]=Lt
Pt=(EstFKVE$Pt)
liPt=Pt-qlils*sqrt(EstFKVEpt$Pt)
lsPt=Pt+qlils*sqrt(EstFKVEpt$Pt)
Ptm=matrix(0,n1,3)
Ptm[,1]=liPt
Ptm[,2]=lsPt
Ptm[,3]=Pt
Wt=(EstFKVE$Wt)
liWt=Wt-qlils*sqrt(EstFKVEpt$Wt)
lsWt=Wt+qlils*sqrt(EstFKVEpt$Wt)
Wtm=matrix(0,n1,3)
Wtm[,1]=liWt
Wtm[,2]=lsWt
Wtm[,3]=Wt

#############################################################################
##
## GRAPH 1:
##
############################################################################  
par(mfrow=c(2,2),mar = c(5.1,4.1,1.3,1.1))
#Et
seq1=seq(1,n1)
seq2=seq(n1,1)
seq3=array(0,c(n1))
require(chron)
Weeks=seq.dates("05/10/11", "04/24/12", by = "weeks") # serie ovos BH
at=seq(as.Date("2011/05/10"), as.Date("2012/04/24"), "weeks")  # Mexer aqui!!
plot(at,Etm[,3],lty=c(2,2,1),ylab="Et",ylim=c(0,max(Etm[2:n1,])),xlab="t (Weeks)",type='l',xaxt = "n")
polygon(c(at, sort(Weeks,decreasing = TRUE)),c((Etm[,1]),rev((Etm[,2]))),col="light grey",border="light grey",angle=90)
par(new=T)
plot(at,Etm[,3],lty=c(2),ylab="Et",xlab="t (Weeks)",ylim=c(0,max(Etm[2:n1,])),type='l',xaxt = "n")
axis.Date(1,at=at,format="%y-%m-%d", col.axis = 1, col = 1)
title("(a)")
#Lt
seq1=seq(1,n1)
seq2=seq(n1,1)
seq3=array(0,c(n1))
plot(at,Ltm[,3],lty=c(2,2,1),ylab="Lt",ylim=c(0,max(Ltm[2:n1,])),xlab="t (Weeks)",type='l',xaxt = "n")
polygon(c(at, sort(Weeks,decreasing = TRUE)),c((Ltm[,1]),rev((Ltm[,2]))),col="light grey",border="light grey",angle=90)
par(new=T)
plot(at,Ltm[,3],lty=c(2),ylab="Lt",ylim=c(0,max(Ltm[2:n1,])),xlab="t (Weeks)",type='l',xaxt = "n")
axis.Date(1,at=at,format="%y-%m-%d", col.axis = 1, col = 1)
title("(b)")
#Pt
seq1=seq(1,n1)
seq2=seq(n1,1)
seq3=array(0,c(n1))
plot(at,Ptm[,3],lty=c(2),ylab="Pt",ylim=c(0,max(Ptm[2:n1,])),xlab="t (Weeks)",type='l',xaxt = "n")
par(new=T)
polygon(c(at, sort(at,decreasing = TRUE)),c((Ptm[,1]),rev((Ptm[,2]))),col="light grey",border="light grey",angle=90)
par(new=T)
plot(at,Ptm[,3],lty=c(2),ylab="Pt",ylim=c(0,max(Ptm[2:n1,])),xlab="t (Weeks)",type='l',xaxt = "n")
axis.Date(1,at=at,format="%y-%m-%d", col.axis = 1, col = 1)
title("(c)")
#Wt
seq1=seq(1,n1)
seq2=seq(n1,1)
seq3=array(0,c(n1))
plot(at,Wtm[,3],lty=c(2),ylab="Wt",ylim=c(0,max(Wtm[2:n1,])),xlab="t (Weeks)",type='l',xaxt = "n")
par(new=T)
polygon(c(at, sort(at,decreasing = TRUE)),c((Wtm[,1]),rev((Wtm[,2]))),col="light grey",border="light grey",angle=90)
par(new=T)
plot(Weeks,Wtm[,3],lty=c(2),ylab="Wt",ylim=c(0,max(Wtm[2:n1,])),xlab="t (Weeks)",type='l',xaxt = "n")
axis.Date(1,at=at,format="%y-%m-%d", col.axis = 1, col = 1)
title("(d)")

#############################################################################
##
## GRAPH 2:
##
############################################################################   
ytt1=EstFKVEytt1$ytt1
Ft=EstFKVEft$Ft
liyt=ytt1-1.65*sqrt(Ft)
lsyt=ytt1+1.65*sqrt(Ft)
liyt
ytt1
lsyt
ytm=matrix(0,n1,4)
liyt[liyt<0]=0
ytm[,1]=liyt
ytm[,2]=lsyt
ytm[,3]=ytt1
ytm[,4]=yt

###############################################################################
windows()
#Et
seq1=seq(1,n1)
seq2=seq(n1,1)
seq3=array(0,c(n1))
require(chron)
Weeks=seq.dates("05/10/11", "04/24/12", by = "weeks")
at=seq(as.Date("2011/05/10"), as.Date("2012/04/24"), "weeks")  
plot(at,ytm[,4],lty=c(1),ylab="yt",ylim=c(0,max(yt)),xlab="t (Weeks)",type='l',xaxt = "n")
par(new=T)
polygon(c(at, sort(at,decreasing = TRUE)),c((ytm[,1]),rev((ytm[,2]))),col="light grey",border="light grey",angle=90)
par(new=T)
plot(at,ytm[,3],lty=c(2),ylab="yt",ylim=c(0,max(yt)),xlab="t (Weeks)",type='l',xaxt = "n")
par(new=T)
plot(at,ytm[,4],lty=c(1),ylab="yt",ylim=c(0,max(yt)),xlab="t (Weeks)",type='l',xaxt = "n")
axis.Date(1,at=at,format="%y-%m-%d", col.axis = 1, col = 1)
#MSE:
mean(abs(yt-Et))
mean(abs(yt-ytt1))

#################################
### Diagnostic
################################
Residuals=(yt-ytt1)
Residuals
Residualsp=(Residuals-mean(Residuals))/sd(Residuals)
Residualsp
Vt=Residualsp
windows()
par(mfrow=c(2,2),mar = c(5.1,4.1,1.3,1.1))
ts.plot(Vt,xlab="t",ylab="Vt",type='l')
title("(a)")
qqnorm(Vt,main="")
title("(b)")
acf(Vt,lag=80,main="",ci=0.99)
title("(c)")
pacf(Vt,lag=80,main="",ci=0.99)
title("(d)")
shapiro.test(Vt)
Box.test(Vt,lag=1)
Box.test(Vt,lag=6)
Box.test(Vt,lag=12)
#############################################################################
