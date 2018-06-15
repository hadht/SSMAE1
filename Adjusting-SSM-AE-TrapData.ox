//////////////////////////////////////////////////////////////////////////////////////////////////////////
//Program for Adjusting a tSate Space Model to Aedes Aegypti Trap Data
//Author: Thiago Rezende
//version: 1.0
//date of creation: 01/01/2018
/////////////////////////////////////////////////////////////////////////////////////////////////////////

//Libraries:
#include <oxstd.h>			 
#include <oxdraw.h>
#import  <maximize>
#import <maxsqp>

//////////////////////////////////////////////////////////////
// declaring (argh!) globals
//////////////////////////////////////////////////////////////
decl g_STSiBoot=500;
decl g_STSiBFGS=100; 
decl g_STSiBurnIn=100;
decl g_STSiMessages=FALSE;
static decl g_STSvYt,ro;
static decl g_STSmPhi,g_ov,g_ro, g_STSmOmega, g_STSmSigma, g_STSmKF;

//	Functions:
 FKMNL(const Yt,const Xt,const sigma2eps,const sig1,const sig2,const sig3,const sig5,const ov,const ind);
 like(const Yt,const Xt,const sigma2eps,const sig1,const sig2,const sig3,const sig5,const ov,const ind);
 STSLoglike2(const vP,const vFunc, const vScore, const mHess);
 STSAdjust2(const vYt,const ind2);
 SimulateMC(const iseed,const size,const tipro);
 STSBootFT(const vYt,const Xt,const sigma2eta,const sigma2eps,const sigma2csi,const ro,const vVt,const vFt,const vKt);
 FuncPrevBoot(const STSvYt,const vYtboot,const sigma2eta,const sigma2csi,const sigma2eps,const ro,const lag);
 STSSimulateNormTransf(const dSigma2Eta,const dSigma2Eps,const dSigma2Csi,const iYtSize,const pro,const Xt,const iseed);
 static decl g_beta,g_Xt;

//////////////////////////////////////////////////////////////////////////////////////
  //KALMAN'S FILTER:
/////////////////////////////////////////////////////////////////////////////////////

FKMNL(const Yt,const Xt,const sigma2eps,const sig1,const sig2,const sig3,const sig5,const ov,const ind){
decl K=10,i,at,att1,ytt1,atoutput,ptoutput,pt,ptt1,vt,Ft,Kt,vtprev,result,att1output,ptt1output,
malfa,valfa,n,h,nPar,vP,zt,Tt,Ttaux,Rt,Qt,ct,Rm,m1,m2,m3,mt4=(1/12.5),ovt=1;

 n=rows(Yt);
 if(n==1){n=rows(Yt');};
 at=zeros(4,1);
 atoutput=zeros(4,n+1);
 att1output=zeros(4,n);
 ptt1output=zeros(4,n);
 ptoutput=zeros(4,n+1);
 att1=zeros(4,1);
 zt=zeros(4,1);
 pt=zeros(4,4);
 ptt1=zeros(4,4);
 vP=zeros(1,4);
 Tt=zeros(4,4);
 Ttaux=zeros(4,4);
 Rt=zeros(4,4);
 ct=zeros(4,1);
 Qt=zeros(4,4);
 vt=zeros(1,n);
 ytt1=zeros(1,n);
 Ft=zeros(1,n);
 Kt=zeros(4,n)';
 nPar=1;
//System Matrices: 
 zt[0]=1;
 zt[1]=0;
 zt[2]=0;
 zt[3]=0;

 Rt[0][0]=1;
 Rt[0][1]=0;
 Rt[0][2]=0;
 Rt[0][3]=0;
 Rt[1][0]=0;
 Rt[1][1]=0;
 Rt[1][2]=0;
 Rt[1][3]=0;
 Rt[2][0]=0;
 Rt[2][1]=0;
 Rt[2][2]=0;
 Rt[2][3]=0;
 Rt[3][0]=0;
 Rt[3][1]=0;
 Rt[3][2]=0;
 Rt[3][3]=0;
 
 //Matrices:
 Qt[0][0]=sig5;
 Qt[0][1]=0;
 Qt[1][1]=0;
 Qt[1][0]=0;
 Qt[2][2]=0;
 Qt[0][2]=0;
 Qt[1][2]=0;
 Qt[2][0]=0;
 Qt[2][1]=0; 
 
//
 ct[0]=0;
 ct[1]=0;
 ct[2]=0;
 ct[3]=0;

//Mortality values:
m1=1/60;
m2=1/2;
m3=1/50;
Rm=((0.5*ov)/(sig1+m1))*(sig1/(sig2+m2))*(sig2/(sig2+m3))*(sig3/mt4);

// Initial Point to specify:////////////////////////////////////////////////
at[0]=K*(1-(1/Rm));	   //Et
at[2]=(sig1/(sig2+m2))*at[0];  //Lt
at[1]=(sig2/(sig2+m3))*at[2];  //Pt
at[3]=(sig3/mt4)*at[1];	//Wt

att1[0]=at[0];
att1[1]=at[1];
att1[2]=at[2];
att1[3]=at[3];
atoutput[0][0]=at[0];
atoutput[1][0]=at[1];
atoutput[2][0]=at[2];
atoutput[3][0]=at[3];

pt[0][0]=50;
pt[0][1]=0;
pt[1][0]=0;
pt[0][2]=0;
pt[2][0]=0;
pt[1][2]=0;
pt[2][1]=0;
pt[1][1]=50;
pt[2][2]=50;
pt[3][3]=50;
pt[3][1]=0;
pt[3][2]=0;
pt[3][0]=0;
pt[0][3]=0;
pt[1][3]=0;
pt[2][3]=0;

 Tt[0][0]=1-1/60-sig1-(ov/(2*K))*at[3]; //Et
 Tt[0][1]=0;
 Tt[0][2]=0;
 Tt[0][3]=(ov/2)*(1-(at[0]/K));
 Tt[1][0]=0;		  //Pt
 Tt[1][1]=1-1/50-sig3;	   
 Tt[1][2]=sig2;
 Tt[1][3]=0;
 Tt[2][0]=sig1;		   //Lt
 Tt[2][1]=0;	   
 Tt[2][2]=1-1/2-sig2;
 Tt[2][3]=0;
 Tt[3][0]=0;			//Wt
 Tt[3][1]=sig3;	   
 Tt[3][2]=0;
 Tt[3][3]=1-mt4;

  //h_t = c_t:	   // nonlinearity:
 ct[0]=(ov/(2*K))*at[0]*at[3];
 ct[1]=0;
 ct[2]=0;
 ct[3]=0;

for(h=0;h<(n);h++){
 //T_t = H_t:	 // nonlinearity:
 Tt[0][0]=1-1/60-sig1-(ov/(2*K))*at[3]; //Et
 Tt[0][1]=0;
 Tt[0][2]=0;
 Tt[0][3]=(ov/2)*(1-(at[0]/K));
 Tt[1][0]=0;		  //Pt
 Tt[1][1]=1-1/50-sig3;	   
 Tt[1][2]=sig2;
 Tt[1][3]=0;
 Tt[2][0]=sig1;		   //Lt
 Tt[2][1]=0;	   
 Tt[2][2]=1-1/2-sig2;
 Tt[2][3]=0;
 Tt[3][0]=0;			//Wt
 Tt[3][1]=sig3;	   
 Tt[3][2]=0;
 Tt[3][3]=1-mt4;


  //h_t = c_t:	   // nonlinearity:
 ct[0]=(ov/(2*K))*at[0]*at[3];
 ct[1]=0;
 ct[2]=0;
 ct[3]=0;
	 //One-step ahead predicted distribution of the states:
	  att1=(Tt*at)+ct;
	  att1output[0][h]=att1[0];
      att1output[1][h]=att1[1];
	  att1output[2][h]=att1[2];
	  att1output[3][h]=att1[3];
	  ptt1=Tt*pt*(Tt')+Rt*Qt*(Rt');
	  ptt1output[0][h]=ptt1[0][0];
	  ptt1output[1][h]=ptt1[1][1];
	  ptt1output[2][h]=ptt1[2][2];
	  ptt1output[3][h]=ptt1[3][3];

   //Quantities of Interest:
	  vt[h]=Yt[h]-(zt')*att1;	// prediction error
	  Ft[h]=(zt')*(ptt1)*zt+sigma2eps;//variance of prediction error
	 decl Ktaux= Tt*ptt1*zt*(1/Ft[h]);
	  Kt[h][0]=Ktaux[0];		//Matrix of Kalman gain
	  Kt[h][1]=Ktaux[1];		
	  Kt[h][2]=Ktaux[2];	
	  Kt[h][3]=Ktaux[3];	
	  ytt1[h]=(zt')*att1;

	//update distribution of the states:
	  //at:
	  at=att1+ptt1*zt*vt[h]*(1/Ft[h]);
	   atoutput[0][h+1]=at[0];
       atoutput[1][h+1]=at[1];
	   atoutput[2][h+1]=at[2];
	   atoutput[3][h+1]=at[3];
	   //Pt:
	  pt=ptt1-ptt1*zt*(1/Ft[h])*(zt')*ptt1;
	  ptoutput[0][h+1]=pt[0][0];
	  ptoutput[1][h+1]=pt[1][1];
	  ptoutput[2][h+1]=pt[2][2];
	  ptoutput[3][h+1]=pt[3][3];

  }
  result=zeros(6,n);
  result[0][:]=vt;
  result[1][:]=Ft;
  result[2:][:]=Kt;
  vtprev=(vt[:]);
  if(ind==1){return(vt);};
  if(ind==2){return(Ft);};
  if(ind==3){return(Kt);};
  if(ind==4){return(atoutput);};
  if(ind==5){return(ptoutput);};
  if(ind==6){return(result);};
  if(ind==7){return((ytt1));};
  if(ind==8){return(zt'*ptt1*zt+sigma2eps);};
  if(ind==9){return(vtprev);};
  if(ind==10){return(att1output);};
  if(ind==11){return(ptt1output);};
  if(ind==12){return(Rm);};

}

 //////////////////////////////////////////////////////////////
// implementing
//////////////////////////////////////////////////////////////
// defining Log-likelihood function
//////////////////////////////////////////////////////////////

like(const Yt,const Xt,const sigma2eps,const sig1,const sig2,const sig3,const sig5,const ov,const ind){
 decl k,vt1,Ft1,Kt1,lc,resul,sig22=0,logft=0;
  decl n=rows(Yt);
  resul=FKMNL(Yt,Xt,sigma2eps,sig1,sig2,sig3,sig5,ov,6);
//  println(resul);
  vt1=resul[0][:];
  Ft1=resul[1][:];
  sig22=(vt1.*vt1)*(1.0./Ft1)';
  logft=(sumr(log(Ft1)));
  lc=-0.5*(sig22+logft);	//log-lik
  if(ind==1){lc=-(n/2)*log(2*3.14)-0.5*(sig22+logft);};
return ((lc));
}

 //////////////////////////////////////////////////////////////
// implementing
//////////////////////////////////////////////////////////////
// defining Log-likelihood function to return values
//////////////////////////////////////////////////////////////
STSLoglike2(const vP,const vFunc, const vScore, const mHess) {
	decl ii, dVar;
     	g_STSmOmega[0][0]=((vP[0]));
		g_STSmOmega[1][1]=((vP[1]));
		g_STSmOmega[2][2]=((vP[2]));
	    g_STSmOmega[3][3]=(vP[3]);
		g_STSmOmega[4][4]=(vP[4]);
		g_ov=((vP[5]));
		vFunc[0]=like(g_STSvYt,g_Xt,g_STSmOmega[4][4],g_STSmOmega[0][0],g_STSmOmega[1][1],
		g_STSmOmega[2][2],g_STSmOmega[3][3],g_ov,0);
	return(1);
}
//////////////////////////////////////////////////////////////
// adjusting the model
//////////////////////////////////////////////////////////////
STSAdjust2(const vYt,const ind2) {
//	var declaration
	decl vP, vFunc, dVar;
	decl ii;
//	get time series
	g_STSvYt=vYt;
    g_STSmOmega=new matrix[5][5];
	g_STSmOmega[0][0]=1/10.5;  //sig1
	g_STSmOmega[1][1]=1/11.7; //sig2
	g_STSmOmega[2][2]=1/4.5; //sig3
	g_STSmOmega[3][3]=1;	//sigmaet
	g_STSmOmega[4][4]=1;	 //sigmaeps
	g_ov=1.0;
//	initialize independent variables
		vP=zeros(6,1);
 	    vP[0]=((3*g_STSmOmega[0][0]));	  
		vP[1]=((3*g_STSmOmega[1][1]));
		vP[2]=((3*g_STSmOmega[2][2]));
		vP[3]=(g_STSmOmega[3][3]);
		vP[4]=(g_STSmOmega[4][4]);
 		vP[5]=((1*(g_ov)));
//	log-likelihood maximization
//	MaxControl(g_STSiBFGS,-1);
	STSLoglike2(vP,&vFunc,0,TRUE);
    MaxControl(g_STSiBFGS,0,0.0000001,0.000005);
//	ii=MaxBFGS(STSLoglike2, &vP, &vFunc, 0, TRUE);
	ii=MaxSQP(STSLoglike2, &vP, &vFunc, 0, TRUE,0,0,<0.01;0.01;0.01;0.01;0.01;0.5>,<1;1;1;1000000;1000000;100>);
	if ((ii!=FALSE)&&(g_STSiMessages==TRUE)) {
		println("STSAdjust:ii=", ii, " ", MaxConvergenceMsg(ii));
	}
	println("STSAdjust:", MaxConvergenceMsg(ii)," using numerical derivatives\n");
	decl result=zeros(4,2);
	result[0][0]= g_STSmOmega[0][0];
	result[0][1]= g_STSmOmega[1][1];
	result[1][0]= g_STSmOmega[2][2];
	result[1][1]= (1/12.5);
	result[2][0]= g_STSmOmega[3][3];
	result[2][1]=g_ov;
	result[3][0]= g_STSmOmega[4][4];
	result[3][1]= 0;
	return(result);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Main program of Application to Real Time Series
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SimulateMC(const iseed,const size,const tipro){

decl Yt,Ytex,Ytboot,Xt,sigma2eta,sigma2eps,i,j,STSomegaboot,RoCov,CsiCov,EtaCov,EpsCov,PrevCov,PrevAssCov,vaux,vt, MC=1,boot=1;
decl ro=zeros(1,MC),Ytomit=zeros(1,MC),linfboot=zeros(1,MC),lsupboot=zeros(1,MC);
decl sigeta=zeros(1,MC),matprev=zeros(1,boot),mptprev=zeros(1,boot);
decl sigcsi=zeros(1,MC),mediabootprev=zeros(1,MC),sig1,sig2,sig3,sig4,sig5,sigeps,ov1;
decl medass=zeros(1,MC),linfass=zeros(1,MC),lsupass=zeros(1,MC);
decl rob=zeros(1,boot),infbro=zeros(1,MC),supbro=zeros(1,MC);
decl sigcsib=zeros(1,boot),infbscsi=zeros(1,MC),supbscsi=zeros(1,MC);
decl sigetab=zeros(1,boot),infbseta=zeros(1,MC),supbseta=zeros(1,MC);
decl sigepsb=zeros(1,boot),infbseps=zeros(1,MC),supbseps=zeros(1,MC);
decl bro=zeros(1,MC);
decl bsigcsi=zeros(1,MC),prevytt1=zeros(1,size),Ftprev=zeros(1,size),att1prev=zeros(4,size),ptt1prev=zeros(4,size);
decl bsigeta=zeros(1,MC);
decl bsigeps=zeros(1,MC);
decl brome=zeros(1,MC);
decl bsigcsime=zeros(1,MC);
decl bsigetame=zeros(1,MC);
decl bsigepsme=zeros(1,MC),at;
decl outp="SSMAE-";
g_STSiBFGS=200;


println("**************Adjusting a tSate Space Model to Aedes Aegypti Trap Data*******************");
fopen(outp~tipro~".out", "la");
println("n=",size);
println("iseed=",iseed);

 CsiCov=EtaCov=EpsCov=RoCov=0;
 ranseed(iseed);
 Xt=zeros(1,size);
 g_Xt=Xt;
 PrevCov=PrevAssCov=0;
 
for(i=0;i<size;i++){  //size-1
println("************************Prediction-t=",i+1,"**********************************************");
 ranseed(i+iseed);
 Yt=log(loadmat("ytovos-costarica.xls"))';
 Yt=(Yt[:(i)])';
 decl n1=rows(Yt);
 g_STSvYt=Yt;
g_STSiBFGS=200;
decl STSomega1=STSAdjust2(Yt,1),atemv;
 sig1=STSomega1[0][0];
 sig2=STSomega1[0][1];
 sig3=STSomega1[1][0];
 sig4=STSomega1[1][1];
 sig5=STSomega1[2][0];
 sigeps=STSomega1[3][0];
 ov1=STSomega1[2][1];
 println("MLE of Static Parameters:",STSomega1);

 //--------------------------------------------------------
at=FKMNL(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,4);
decl pt1=FKMNL(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,5);
vt=FKMNL(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,1);
decl Rmm=FKMNL(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,12);
//println("R0=",Rmm);

//like(const Yt,const Xt,const sigma2eps,const sig1,const sig2,const sig3,const sig5,const ov,const ind);
decl likv=like(Yt,g_Xt, sigeps,0.220, 0.085, 0.217,sig5, 1, 0);
//println("likv-flr=",likv);
likv=like(Yt,g_Xt, sigeps,0.303, 0.125,  0.323,sig5, 1, 0);
//println("likv-fir=",likv);
likv=like(Yt,g_Xt, sigeps,0.500, 0.227,  0.167,sig5, 1, 0);
//println("likv-fhr=",likv);

savemat("atEggDengueModelv21.xls", at',{"Et","Pt","lt","Wt"});
savemat("atEggDengueModelv21.csv", at',{"Et","Pt","lt","Wt"});
savemat("ptEggDengueModelv21.xls", pt1',{"Et","Pt","lt","Wt"});
savemat("ptEggDengueModelv21.csv", pt1',{"Et","Pt","lt","Wt"});
savemat("vtEggDengueModelv21.xls", vt',{"vt"});
savemat("vtEggDengueModelv21.csv", vt',{"vt"});

decl ptemv=FKMNL(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,5);
     atemv=FKMNL(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,4);
decl vtemv=FKMNL(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,1);

decl att1emv=FKMNL(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,10);
decl ptt1emv=FKMNL(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,11);

att1prev[0:3][i]=att1emv[0:3][i];
ptt1prev[0:3][i]=ptt1emv[0:3][i];

savemat("att1EggDengueModelv21.xls", att1prev',{"Et","Pt","lt","Wt"});
savemat("att1EggDengueModelv21.csv", att1prev',{"Et","Pt","lt","Wt"});
savemat("ptt1EggDengueModelv21.xls",ptt1prev',{"Et","Pt","lt","Wt"});
savemat("ptt1EggDengueModelv21.csv",ptt1prev',{"Et","Pt","lt","Wt"});

decl ytt1=FKMNL(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,7);
decl Ft=FKMNL(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,2);
savemat("ytt1EggDengueModelv21.xls", ytt1',{"ytt1"});
savemat("ytt1EggDengueModelv21.csv", ytt1',{"ytt1"});
savemat("FtEggDengueModelv21.xls", Ft',{"Ft"});
savemat("FtEggDengueModelv21.csv", Ft',{"Ft"});

prevytt1[i]=ytt1[i];
Ftprev[i]=Ft[i];

} //fim MC

println("**************************Forecast**************************************************");
println("meanprevass=",meanr(medass));
println("Meanlinfass=",meanr(linfass));
println("Meansupass=",meanr(lsupass));
println("liminfass=",linfass);
println("limsupass=",lsupass);
println("Coverage=",PrevAssCov/MC);



println("*****************Results MC****************************************");
println("Oviposicao:",ov1);
println("sig1:(modelo csi)",sig1);
println("sig2:(modelo csi)",sig2);
println("sig3:(modelo csi)",sig3);
println("mt_4:(modelo csi)",sig4);
println("sigma2eta:",sig5);
println("sigma2eps:",sigeps);
println("Razao-Sinal-Ruido:",sig5/sigeps);

println("-----------------------------------------------------");
   println("Adjust Statistics:");
   decl vFunc=like(Yt,g_Xt,(sigeps),(sig1),(sig2),(sig3),(sig5),ov1,1);
   decl size1=size-1;
   println("size=",size1);
   println("vFunc=",vFunc);
   println("-2*logvero=",-2*vFunc[0]);
   println("STSAdjust:AIC=",(1/size1)*(-2*vFunc[0]+2*(3+3)));

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Main Program
///////////////////////////////////////////////////////////////////////////////////////////////////////////

main(){
   println("----------------------------------------------");
   decl time=timer();
   SimulateMC(65,51,"AdjustingSSMEggTimeSeries");
   println("Tempo Computacional Geral:");
   println(timespan(time,timer()));

}