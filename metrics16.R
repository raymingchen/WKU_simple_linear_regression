######### Metrics based on PCA
rm(list=ls());
library("scatterplot3d");
##Delta 1 metric
##B1,B2 are n-by-n matrices whose rows are basis vectors
Delta1=function(p_colMat,B1,q_colMat,B2)  
{
n=dim(p_colMat)[1];
difVec=B1%*%p_colMat-B2%*%q_colMat;
d1=t(difVec)%*%(difVec);

d2=0;
for(i in 1:n)
  { 
  d2=d2+(B1[i,]-B2[i,])%*%(B1[i,]-B2[i,]);
  }
delta1=sqrt(d1)+1/n*sqrt(d2);
return(delta1);
}

##Dmat is a m1-by-n data matrix, Emat is a m2-by-n data matrix
##where m1,m2 are number of data points; n is the dimension
##B1,B2 are n-by-n matrices whose rows are basic vectors
DisMat=function(Dmat,Emat,B1,B2) 
{
m1=dim(Dmat)[1];
m2=dim(Emat)[2];
disMat=matrix(nrow=m1,ncol=m2);
for(i in 1:m1)
{
p_colMat=t(Dmat[i,,drop=FALSE]);
    for(j in 1:m2)
    {
    q_colMat=t(Emat[j,,drop=FALSE]);
    delta1=Delta1(p_colMat,B1,q_colMat,B2); 
    disMat[i,j]=delta1;
    }
}
return(disMat);
}

##Hausdorff metric
Hausdorff=function(disMat) ##disMat: a m-by-n distance matrix
{
min_row=apply(disMat,1,min);
min_col=apply(disMat,2,min);
avg_row=mean(min_row);
avg_col=mean(min_col);
hausdorff=max(avg_row,avg_col);
return(hausdorff);
}

##Delta2 metric
Delta2=function(Dmat,Emat,B1,B2) 
{
disMat=DisMat(Dmat,Emat,B1,B2);
hausdorff=Hausdorff(disMat);
return(hausdorff);
}


Pca_basis=function(Dmat)
{
eig=eigen(t(Dmat)%*%Dmat);
eigVec=eig$vectors;
pca_basis=t(eigVec);  ## after the transpose, the row vectors are the basic vector
return(pca_basis);
}

##Delta3 metric
##alpha,beta are positive reals for weighing
Delta3=function(p_colMat,Dmat,q_colMat,Emat,alpha,beta) 
{
B1=Pca_basis(Dmat);
B2=Pca_basis(Emat);
delta1=Delta1(p_colMat,B1,q_colMat,B2);
delta2=Delta2(Dmat,Emat,B1,B2);
delta3=alpha*delta1+beta*delta2;
return(delta3);
}


##Delta0 Mahalanobis metric
##clVec:column vector
##DMat is a m(data points)-by-n(dimension) matrix
MahaMetric=function(clVec,DMat)  
{
X=clVec;
X_demean=clVec-mean(X);
pca_basis=Pca_basis(Dmat);
mahaMetric=(t(X_demean)%*%solve(pca_basis)%*%X_demean);
return(mahaMetric);
}

##--compute the distance between featured matrices (row-wise distances)
FeatureDis=function(mat1,mat2)  ##mat1,2:m-by-n matrices
{
m=dim(mat1)[1];
n=dim(mat2)[2];
dif=mat1-mat2;
featureDis=sum(sqrt(apply((dif*dif),1,sum)));
return(featureDis);
}


##---Mahalanobis dis1
Mahalanobis=function(DD,num)  ##m-by-n matrix; num:index; n:dimension, m:#data points
{
n=dim(DD)[2];
DD_inc=DD;
DD_exc=DD[-c(num),];  ##delete num'th row of DD 
DDmean_inc=apply(DD_inc,2,mean);
DDmean_exc=apply(DD_exc,2,mean);
Vec=DD[num,];   ## the outlined vector
Vec_demean=Vec-DDmean_inc;
eigVa_inc=eigen(t(DD_inc)%*%DD_inc)$values;
eigVa_exc=eigen(t(DD_exc)%*%DD_exc)$values;
invS_inc=1/eigVa_inc*diag(n);
invS_exc=1/eigVa_exc*diag(n);
invCov_inc=solve(cov(DD_inc));
invCov_exc=solve(cov(DD_exc));
invCov_inc=solve(cov(DD_inc));
manhalanobis_inc=sqrt((t(Vec_demean)%*%invCov_inc)%*%Vec_demean);
manhalanobis_exc=sqrt((t(Vec_demean)%*%invCov_exc)%*%Vec_demean);
return(c(manhalanobis_inc,manhalanobis_exc));
}


##--Dependcy and Integration
DepInt=function(DD,num)  ##DD:m-by-n matrix
{
vector=DD[num,];
DD_inc=DD;
DD_exc=DD[-c(num),];  ##delete num'th row of DD 
eigVa_inc=eigen(t(DD_inc)%*%DD_inc)$values;
eigVa_exc=eigen(t(DD_exc)%*%DD_exc)$values;
eigVec_inc=eigen(t(DD_inc)%*%DD_inc)$vectors;
eigVec_exc=eigen(t(DD_exc)%*%DD_exc)$vectors;
norm_vector=sqrt(vector%*%vector);
cos_inc=(eigVec_inc%*%vector)/as.numeric(norm_vector);
cos_exc=(eigVec_exc%*%vector)/as.numeric(norm_vector);
weight_inc=eigVa_inc/sum(eigVa_inc);
weight_exc=eigVa_exc/sum(eigVa_exc);
dependency=weight_inc%*%cos_inc;
integration=weight_exc%*%cos_exc;
return(c(dependency,integration));
}


###===========================Data Analysis=====================================================%%

data1_MAIN_TW_HK=read.csv(file.choose());
data2_MAIN_TW_HK=read.csv(file.choose());
data3_MAIN_TW_HK=read.csv(file.choose());

data11_MAIN_TW_HK=data1_MAIN_TW_HK[1:33,c(1,8:25)];
data22_MAIN_TW_HK=data2_MAIN_TW_HK[1:33,c(1,8:25)];
data33_MAIN_TW_HK=data3_MAIN_TW_HK[1:33,c(1,8:25)];

data111_MAIN_TW_HK=data11_MAIN_TW_HK[1:33,c(2:19)];
data222_MAIN_TW_HK=data22_MAIN_TW_HK[1:33,c(2:19)];
data333_MAIN_TW_HK=data33_MAIN_TW_HK[1:33,c(2:19)];
data111_MAIN_TW=data11_MAIN_TW_HK[1:32,c(2:19)];
data222_MAIN_TW=data22_MAIN_TW_HK[1:32,c(2:19)];
data333_MAIN_TW=data33_MAIN_TW_HK[1:32,c(2:19)];
data111_MAIN_HK=data11_MAIN_TW_HK[c(1:31,33),c(2:19)];
data222_MAIN_HK=data22_MAIN_TW_HK[c(1:31,33),c(2:19)];
data333_MAIN_HK=data33_MAIN_TW_HK[c(1:31,33),c(2:19)];
data111_MAIN=data11_MAIN_TW_HK[1:31,c(2:19)];
data222_MAIN=data22_MAIN_TW_HK[1:31,c(2:19)];
data333_MAIN=data33_MAIN_TW_HK[1:31,c(2:19)];


dataPoints_MAIN_TW_HK=dim(data111_MAIN_TW_HK)[1];
dataPoints_MAIN_TW=dim(data111_MAIN_TW)[1];
dataPoints_MAIN_HK=dim(data111_MAIN_HK)[1];
dataPoints_MAIN=dim(data111_MAIN)[1];

times=dim(data111_MAIN_TW_HK)[2];

dataList_MAIN_TW_HK=as.list(numeric(times));dim(dataList_MAIN_TW_HK)=c(1,times);
dataList_MAIN_TW=as.list(numeric(times));dim(dataList_MAIN_TW)=c(1,times);
dataList_MAIN_HK=as.list(numeric(times));dim(dataList_MAIN_HK)=c(1,times);
dataList_MAIN=as.list(numeric(times));dim(dataList_MAIN)=c(1,times);

PcaList_MAIN_TW_HK=as.list(numeric(times));dim(PcaList_MAIN_TW_HK)=c(1,times);
PcaList_MAIN_TW=as.list(numeric(times));dim(PcaList_MAIN_TW)=c(1,times);
PcaList_MAIN_HK=as.list(numeric(times));dim(PcaList_MAIN_HK)=c(1,times);
PcaList_MAIN=as.list(numeric(times));dim(PcaList_MAIN)=c(1,times);

for(i in 1:times)
{
mat1=matrix(nrow=dataPoints_MAIN_TW_HK,ncol=3);  ##3 industries
mat2=matrix(nrow=dataPoints_MAIN_TW,ncol=3);  ##3 industries
mat3=matrix(nrow=dataPoints_MAIN_HK,ncol=3);  ##3 industries
mat4=matrix(nrow=dataPoints_MAIN,ncol=3);  ##3 industries

mat1[,1]=data111_MAIN_TW_HK[,i];
mat1[,2]=data222_MAIN_TW_HK[,i];
mat1[,3]=data333_MAIN_TW_HK[,i];
mat2[,1]=data111_MAIN_TW[,i];
mat2[,2]=data222_MAIN_TW[,i];
mat2[,3]=data333_MAIN_TW[,i];
mat3[,1]=data111_MAIN_HK[,i];
mat3[,2]=data222_MAIN_HK[,i];
mat3[,3]=data333_MAIN_HK[,i];
mat4[,1]=data111_MAIN[,i];
mat4[,2]=data222_MAIN[,i];
mat4[,3]=data333_MAIN[,i];

dataList_MAIN_TW_HK[[1,i]]=mat1;  ## the row vectors are the basic vector
dataList_MAIN_TW[[1,i]]=mat2;  ## the row vectors are the basic vector
dataList_MAIN_HK[[1,i]]=mat3;  ## the row vectors are the basic vector
dataList_MAIN[[1,i]]=mat4;  ## the row vectors are the basic vector

PcaList_MAIN_TW_HK[[i]]=Pca_basis(mat1);
PcaList_MAIN_TW[[i]]=Pca_basis(mat2);
PcaList_MAIN_HK[[i]]=Pca_basis(mat3);
PcaList_MAIN[[i]]=Pca_basis(mat4);
}


len_MAIN_TW_HK=length(PcaList_MAIN_TW_HK)-1;
len_MAIN_TW=length(PcaList_MAIN_TW)-1;
len_MAIN_HK=length(PcaList_MAIN_HK)-1;
len_MAIN=length(PcaList_MAIN)-1;

vec_MAIN_TW_HK=rep(0,len_MAIN_TW_HK);
vec_MAIN_TW=rep(0,len_MAIN_TW);
vec_MAIN_HK=rep(0,len_MAIN_HK);
vec_MAIN=rep(0,len_MAIN);

for(t in 1:len_MAIN_TW_HK)
{
featureDis_MAIN_TW_HK=FeatureDis(PcaList_MAIN_TW_HK[[t]],PcaList_MAIN_TW_HK[[t+1]]); 
vec_MAIN_TW_HK[t]=featureDis_MAIN_TW_HK;
}
for(t in 1:len_MAIN_TW)
{
featureDis_MAIN_TW=FeatureDis(PcaList_MAIN_TW[[t]],PcaList_MAIN_TW_HK[[t+1]]); 
vec_MAIN_TW[t]=featureDis_MAIN_TW;
}
for(t in 1:len_MAIN_HK)
{
featureDis_MAIN_HK=FeatureDis(PcaList_MAIN_HK[[t]],PcaList_MAIN_TW_HK[[t+1]]); 
vec_MAIN_HK[t]=featureDis_MAIN_HK;
}
for(t in 1:len_MAIN)
{
featureDis_MAIN=FeatureDis(PcaList_MAIN[[t]],PcaList_MAIN[[t+1]]); 
vec_MAIN[t]=featureDis_MAIN;
}


dependencyMat=matrix(nrow=33,ncol=18);
integrationMat=matrix(nrow=33,ncol=18);

for(t in 1:18)
{
DD=dataList_MAIN_TW_HK[[1,t]];
   for(p in 1:33)
   {
   num=p;
   depInt=DepInt(DD,num);
   dependencyMat[p,t]=depInt[1];
   integrationMat[p,t]=depInt[2];
   }
}

dependency_TW=dependencyMat[32,];
dependency_HK=dependencyMat[33,];
integration_TW=integrationMat[32,];
integration_HK=integrationMat[33,];

apply(-dependencyMat,2,rank);
apply(-integrationMat,2,rank);

mahalanobis_inc_TW=DepInt(DD,32)[1];
mahalanobis_exc_TW=DepInt(DD,32)[2];
mahalanobis_inc_HK=DepInt(DD,33)[1];
mahalanobis_exc_HK=DepInt(DD,33)[2];

wt_inc_TW=mahalanobis_inc_TW/(mahalanobis_inc_TW+mahalanobis_inc_HK);
wt_inc_HK=mahalanobis_inc_HK/(mahalanobis_inc_TW+mahalanobis_inc_HK);
wt_exc_TW=mahalanobis_exc_TW/(mahalanobis_exc_TW+mahalanobis_exc_HK);
wt_exc_HK=mahalanobis_exc_HK/(mahalanobis_exc_TW+mahalanobis_exc_HK);

dependency_TWHK=wt_exc_TW*dependency_TW+wt_exc_HK*dependency_HK;
integration_TWHK=wt_inc_TW*integration_TW+wt_inc_HK*integration_HK;

par(mfrow=c(3,2));
scale=2003:2020;
plot(scale,dependency_TW,xlab="Year",ylab="Dependency",main="Dependency for Taiwan");
plot(scale,integration_TW,xlab="Year",ylab="Integration",main="Integration for Taiwan");
plot(scale,dependency_HK,xlab="Year",ylab="Dependency",main="Dependency for Hong Kong");
plot(scale,integration_HK,xlab="Year",ylab="Integration",main="Integration for Hong Kong");
plot(scale,dependency_TWHK,xlab="Year",ylab="Dependency",main="Dependency for Taiwan and Hong Kong");
plot(scale,integration_TWHK,xlab="Year",ylab="Integration",main="Integration for Taiwan and Hong Kong");



##---regional dependency and integration

TER=data1_MAIN_TW_HK[1:31,3];
Eastern=which(TER=="Eastern Regions");
Central=which(TER=="Central Regions");
Western=which(TER=="Western Regions");
Eastern_TW=union(Eastern,32);
Eastern_HK=union(Eastern,33);
Eastern_TW_HK=union(union(Eastern,32),33);
Central_TW=union(Central,32);
Central_HK=union(Central,33);
Central_TW_HK=union(union(Central,32),33);
Western_TW=union(Western,32);
Western_HK=union(Western,33);
Western_TW_HK=union(union(Western,32),33);

Dep_TW_Eastern=rep(0,18);
Int_TW_Eastern=rep(0,18);
Dep_HK_Eastern=rep(0,18);
Int_HK_Eastern=rep(0,18);
for(t in 1:18)
{
dataList_MAIN_TW_Eastern=dataList_MAIN_TW_HK[[t]];
dataList_MAIN_HK_Eastern=dataList_MAIN_TW_HK[[t]];
DD_TW_Eastern=dataList_MAIN_TW_Eastern[Eastern_TW,];
DD_HK_Eastern=dataList_MAIN_HK_Eastern[Eastern_HK,];
Dep_TW_Eastern[t]=DepInt(DD_TW_Eastern,length(Eastern_TW))[1]; 
Int_TW_Eastern[t]=DepInt(DD_TW_Eastern,length(Eastern_TW))[2]; 
Dep_HK_Eastern[t]=DepInt(DD_HK_Eastern,length(Eastern_HK))[1]; 
Int_HK_Eastern[t]=DepInt(DD_HK_Eastern,length(Eastern_HK))[2]; 
}

Dep_TW_Central=rep(0,18);
Int_TW_Central=rep(0,18);
Dep_HK_Central=rep(0,18);
Int_HK_Central=rep(0,18);
for(t in 1:18)
{
dataList_MAIN_TW_Central=dataList_MAIN_TW_HK[[t]];
dataList_MAIN_HK_Central=dataList_MAIN_TW_HK[[t]];
DD_TW_Central=dataList_MAIN_TW_Central[Central_TW,];
DD_HK_Central=dataList_MAIN_HK_Central[Central_HK,];
Dep_TW_Central[t]=DepInt(DD_TW_Central,length(Central_TW))[1]; 
Int_TW_Central[t]=DepInt(DD_TW_Central,length(Central_TW))[2]; 
Dep_HK_Central[t]=DepInt(DD_HK_Central,length(Central_HK))[1]; 
Int_HK_Central[t]=DepInt(DD_HK_Central,length(Central_HK))[2]; 
}

Dep_TW_Western=rep(0,18);
Int_TW_Western=rep(0,18);
Dep_HK_Western=rep(0,18);
Int_HK_Western=rep(0,18);
for(t in 1:18)
{
dataList_MAIN_TW_Western=dataList_MAIN_TW_HK[[t]];
dataList_MAIN_HK_Western=dataList_MAIN_TW_HK[[t]];
DD_TW_Western=dataList_MAIN_TW_Western[Western_TW,];
DD_HK_Western=dataList_MAIN_HK_Western[Western_HK,];
Dep_TW_Western[t]=DepInt(DD_TW_Western,length(Western_TW))[1]; 
Int_TW_Western[t]=DepInt(DD_TW_Western,length(Western_TW))[2]; 
Dep_HK_Western[t]=DepInt(DD_HK_Western,length(Western_HK))[1]; 
Int_HK_Western[t]=DepInt(DD_HK_Western,length(Western_HK))[2]; 
}

par(mfrow=c(2,2))
scale=2003:2020;
plot(scale,Dep_TW_Eastern,    lwd=2,type="b",xlab="Year",ylab="Dependency",col="red"  ,main="Regional Dependency for TW",ylim=c(0.15,1));
lines(scale,Dep_TW_Central,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="blue" ,main="Regional Dependency for TW");
lines(scale,Dep_TW_Western,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="green",main="Regional Dependency for TW");
plot(scale,Dep_HK_Eastern,    lwd=2,type="b",xlab="Year",ylab="Dependency",col="red"  ,main="Regional Dependency for HK",ylim=c(0.15,1));
lines(scale,Dep_HK_Central,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="blue" ,main="Regional Dependency for HK");
lines(scale,Dep_HK_Western,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="green",main="Regional Dependency for HK");
plot(scale,Int_TW_Eastern,    lwd=2,type="b",xlab="Year",ylab="Dependency",col="red"  ,main="Regional Integration for TW",ylim=c(-0.7,1));
lines(scale,Int_TW_Central,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="blue" ,main="Regional Integration for TW");
lines(scale,Int_TW_Western,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="green",main="Regional Integration for TW");
plot(scale,Int_HK_Eastern,    lwd=2,type="b",xlab="Year",ylab="Dependency",col="red"  ,main="Regional Integration for HK",ylim=c(-0.7,1));
lines(scale,Int_HK_Central,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="blue" ,main="Regional Integration for HK");
lines(scale,Int_HK_Western,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="green",main="Regional Integration for HK");

##---general dependency and integration
General=data1_MAIN_TW_HK[1:31,4];
Huabei=which(General=="Huabei");
Dongbei=which(General=="Dongbei");
Huadong=which(General=="Huadong");
Zhongnan=which(General=="Zhongnan");
Xinan=which(General=="Xinan");
Xibei=which(General=="Xibei");

Huabei_TW=union(Huabei,32);
Huabei_HK=union(Huabei,33);
Dongbei_TW=union(Dongbei,32);
Dongbei_HK=union(Dongbei,33);
Huadong_TW=union(Huadong,32);
Huadong_HK=union(Huadong,33);
Zhongnan_TW=union(Zhongnan,32);
Zhongnan_HK=union(Zhongnan,33);
Xinan_TW=union(Xinan,32);
Xinan_HK=union(Xinan,33);
Xibei_TW=union(Xibei,32);
Xibei_HK=union(Xibei,33);

Dep_TW_Huabei=rep(0,18);
Int_TW_Huabei=rep(0,18);
Dep_TW_Dongbei=rep(0,18);
Int_TW_Dongbei=rep(0,18);
Dep_TW_Huadong=rep(0,18);
Int_TW_Huadong=rep(0,18);
Dep_TW_Zhongnan=rep(0,18);
Int_TW_Zhongnan=rep(0,18);
Dep_TW_Xinan=rep(0,18);
Int_TW_Xinan=rep(0,18);
Dep_TW_Xibei=rep(0,18);
Int_TW_Xibei=rep(0,18);

Dep_HK_Huabei=rep(0,18);
Int_HK_Huabei=rep(0,18);
Dep_HK_Dongbei=rep(0,18);
Int_HK_Dongbei=rep(0,18);
Dep_HK_Huadong=rep(0,18);
Int_HK_Huadong=rep(0,18);
Dep_HK_Zhongnan=rep(0,18);
Int_HK_Zhongnan=rep(0,18);
Dep_HK_Xinan=rep(0,18);
Int_HK_Xinan=rep(0,18);
Dep_HK_Xibei=rep(0,18);
Int_HK_Xibei=rep(0,18);


for(t in 1:18)
{
dataList_MAIN_TW_Huabei=dataList_MAIN_TW_HK[[t]];
dataList_MAIN_HK_Huabei=dataList_MAIN_TW_HK[[t]];

DD_TW_Huabei=dataList_MAIN_TW_Huabei[Huabei_TW,];
DD_HK_Huabei=dataList_MAIN_HK_Huabei[Huabei_HK,];
DD_TW_Dongbei=dataList_MAIN_TW_Huabei[Dongbei_TW,];
DD_HK_Dongbei=dataList_MAIN_HK_Huabei[Dongbei_HK,];
DD_TW_Huadong=dataList_MAIN_TW_Huabei[Huadong_TW,];
DD_HK_Huadong=dataList_MAIN_HK_Huabei[Huadong_HK,];
DD_TW_Zhongnan=dataList_MAIN_TW_Huabei[Zhongnan_TW,];
DD_HK_Zhongnan=dataList_MAIN_HK_Huabei[Zhongnan_HK,];
DD_TW_Xinan=dataList_MAIN_TW_Huabei[Xinan_TW,];
DD_HK_Xinan=dataList_MAIN_HK_Huabei[Xinan_HK,];
DD_TW_Xibei=dataList_MAIN_TW_Huabei[Xibei_TW,];
DD_HK_Xibei=dataList_MAIN_HK_Eastern[Xibei_HK,];

Dep_TW_Huabei[t]=DepInt(DD_TW_Huabei,dim(DD_TW_Huabei)[1])[1]; 
Int_TW_Huabei[t]=DepInt(DD_TW_Huabei,dim(DD_TW_Huabei)[1])[2]; 
Dep_TW_Dongbei[t]=DepInt(DD_TW_Dongbei,dim(DD_TW_Dongbei)[1])[1]; 
Int_TW_Dongbei[t]=DepInt(DD_TW_Dongbei,dim(DD_TW_Dongbei)[1])[2]; 
Dep_TW_Huadong[t]=DepInt(DD_TW_Huadong,dim(DD_TW_Huadong)[1])[1]; 
Int_TW_Huadong[t]=DepInt(DD_TW_Huadong,dim(DD_TW_Huadong)[1])[2]; 
Dep_TW_Zhongnan[t]=DepInt(DD_TW_Zhongnan,dim(DD_TW_Zhongnan)[1])[1]; 
Int_TW_Zhongnan[t]=DepInt(DD_TW_Zhongnan,dim(DD_TW_Zhongnan)[1])[2]; 
Dep_TW_Xinan[t]=DepInt(DD_TW_Xinan,dim(DD_TW_Xinan)[1])[1]; 
Int_TW_Xinan[t]=DepInt(DD_TW_Xinan,dim(DD_TW_Xinan)[1])[2]; 
Dep_TW_Xibei[t]=DepInt(DD_TW_Xibei,dim(DD_TW_Xibei)[1])[1]; 
Int_TW_Xibei[t]=DepInt(DD_TW_Xibei,dim(DD_TW_Xibei)[1])[2]; 

Dep_HK_Huabei[t]=DepInt(DD_HK_Huabei,dim(DD_HK_Huabei)[1])[1]; 
Int_HK_Huabei[t]=DepInt(DD_HK_Huabei,dim(DD_HK_Huabei)[1])[2]; 
Dep_HK_Dongbei[t]=DepInt(DD_HK_Dongbei,dim(DD_HK_Dongbei)[1])[1]; 
Int_HK_Dongbei[t]=DepInt(DD_HK_Dongbei,dim(DD_HK_Dongbei)[1])[2]; 
Dep_HK_Huadong[t]=DepInt(DD_HK_Huadong,dim(DD_HK_Huadong)[1])[1]; 
Int_HK_Huadong[t]=DepInt(DD_HK_Huadong,dim(DD_HK_Huadong)[1])[2]; 
Dep_HK_Zhongnan[t]=DepInt(DD_HK_Zhongnan,dim(DD_HK_Zhongnan)[1])[1]; 
Int_HK_Zhongnan[t]=DepInt(DD_HK_Zhongnan,dim(DD_HK_Zhongnan)[1])[2]; 
Dep_HK_Xinan[t]=DepInt(DD_HK_Xinan,dim(DD_HK_Xinan)[1])[1]; 
Int_HK_Xinan[t]=DepInt(DD_HK_Xinan,dim(DD_HK_Xinan)[1])[2]; 
Dep_HK_Xibei[t]=DepInt(DD_HK_Xibei,dim(DD_HK_Xibei)[1])[1]; 
Int_HK_Xibei[t]=DepInt(DD_HK_Xibei,dim(DD_HK_Xibei)[1])[2]; 
}


par(mfrow=c(2,2))
scale=2003:2020;
plot(scale,Dep_TW_Huabei,    lwd=2,type="b",xlab="Year",ylab="Dependency",col="red"  ,main="General Dependency for TW",ylim=c(0,1));
lines(scale,Dep_TW_Dongbei,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="blue" ,main="General Dependency for TW");
lines(scale,Dep_TW_Huadong,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="green",main="General Dependency for TW");
lines(scale,Dep_TW_Zhongnan,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="yellow" ,main="General Dependency for TW");
lines(scale,Dep_TW_Xinan,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="violet",main="General Dependency for TW");
lines(scale,Dep_TW_Xibei,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="black" ,main="General Dependency for TW");

plot(scale,Dep_HK_Huabei,    lwd=2,type="b",xlab="Year",ylab="Dependency",col="red"  ,main="General Dependency for HK",ylim=c(0.65,1));
lines(scale,Dep_HK_Dongbei,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="blue" ,main="General Dependency for HK");
lines(scale,Dep_HK_Huadong,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="green",main="General Dependency for HK");
lines(scale,Dep_HK_Zhongnan,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="yellow" ,main="General Dependency for HK");
lines(scale,Dep_HK_Xinan,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="violet",main="General Dependency for HK");
lines(scale,Dep_HK_Xibei,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="black" ,main="General Dependency for HK");

plot(scale,Int_TW_Huabei,    lwd=2,type="b",xlab="Year",ylab="Integration",col="red"  ,main="General Integration for TW",ylim=c(-0.5,1));
lines(scale,Int_TW_Dongbei,   lwd=2,type="b",xlab="Year",ylab="Integration",col="blue" ,main="General Integration for TW");
lines(scale,Int_TW_Huadong,   lwd=2,type="b",xlab="Year",ylab="Integration",col="green",main="General Integration for TW");
lines(scale,Int_TW_Zhongnan,   lwd=2,type="b",xlab="Year",ylab="Integration",col="yellow" ,main="General Integration for TW");
lines(scale,Int_TW_Xinan,   lwd=2,type="b",xlab="Year",ylab="Integration",col="violet",main="General Integration for TW");
lines(scale,Int_TW_Xibei,   lwd=2,type="b",xlab="Year",ylab="Integration",col="black" ,main="General Integration for TW");

plot(scale,Int_HK_Huabei,    lwd=2,type="b",xlab="Year",ylab="Integration",col="red"  ,main="General Integration for HK",ylim=c(-0.5,1));
lines(scale,Int_HK_Dongbei,   lwd=2,type="b",xlab="Year",ylab="Integration",col="blue" ,main="General Integration for HK");
lines(scale,Int_HK_Huadong,   lwd=2,type="b",xlab="Year",ylab="Integration",col="green",main="General Integration for HK");
lines(scale,Int_HK_Zhongnan,   lwd=2,type="b",xlab="Year",ylab="Integration",col="yellow" ,main="General Integration for HK");
lines(scale,Int_HK_Xinan,   lwd=2,type="b",xlab="Year",ylab="Integration",col="violet",main="General Integration for HK");
lines(scale,Int_HK_Xibei,   lwd=2,type="b",xlab="Year",ylab="Integration",col="black" ,main="General Integration for HK");

##---EER dependency and integration
EER=data1_MAIN_TW_HK[1:31,6];
NorthernCoastalAreas=which(EER=="Northern Coastal Areas");
MiddleReachesoftheYellowRiver=which(EER=="Middle Reaches of the Yellow River");
Northeast=which(EER=="Northeast ");
EasternCoastalAreas=which(EER=="Eastern Coastal Areas");
MiddleReachesoftheYangzseRiver=which(EER=="Middle Reaches of the Yangzse River");
SouthernCoastalAreas=which(EER=="Southern Coastal Areas");
SoutheastZone=which(EER=="Southeast Zone");
NorthwestZone=which(EER=="Northwest Zone");

NorthernCoastalAreas_TW=union(NorthernCoastalAreas,32);
NorthernCoastalAreas_HK=union(NorthernCoastalAreas,33);
MiddleReachesoftheYellowRiver_TW=union(MiddleReachesoftheYellowRiver,32);
MiddleReachesoftheYellowRiver_HK=union(MiddleReachesoftheYellowRiver,33);
Northeast_TW=union(Northeast,32);
Northeast_HK=union(Northeast,33);
EasternCoastalAreas_TW=union(EasternCoastalAreas,32);
EasternCoastalAreas_HK=union(EasternCoastalAreas,33);
MiddleReachesoftheYangzseRiver_TW=union(MiddleReachesoftheYangzseRiver,32);
MiddleReachesoftheYangzseRiver_HK=union(MiddleReachesoftheYangzseRiver,33);
SouthernCoastalAreas_TW=union(SouthernCoastalAreas,32);
SouthernCoastalAreas_HK=union(SouthernCoastalAreas,33);
SoutheastZone_TW=union(SoutheastZone,32);
SoutheastZone_HK=union(SoutheastZone,33);
NorthwestZone_TW=union(NorthwestZone,32);
NorthwestZone_HK=union(NorthwestZone,33);

Dep_TW_NorthernCoastalAreas=rep(0,18);
Int_TW_NorthernCoastalAreas=rep(0,18);
Dep_TW_MiddleReachesoftheYellowRiver=rep(0,18);
Int_TW_MiddleReachesoftheYellowRiver=rep(0,18);
Dep_TW_Northeast=rep(0,18);
Int_TW_Northeast=rep(0,18);
Dep_TW_EasternCoastalAreas=rep(0,18);
Int_TW_EasternCoastalAreas=rep(0,18);
Dep_TW_MiddleReachesoftheYangzseRiver=rep(0,18);
Int_TW_MiddleReachesoftheYangzseRiver=rep(0,18);
Dep_TW_SouthernCoastalAreas=rep(0,18);
Int_TW_SouthernCoastalAreas=rep(0,18);
Dep_TW_SoutheastZone=rep(0,18);
Int_TW_SoutheastZone=rep(0,18);
Dep_TW_NorthwestZone=rep(0,18);
Int_TW_NorthwestZone=rep(0,18);

Dep_HK_NorthernCoastalAreas=rep(0,18);
Int_HK_NorthernCoastalAreas=rep(0,18);
Dep_HK_MiddleReachesoftheYellowRiver=rep(0,18);
Int_HK_MiddleReachesoftheYellowRiver=rep(0,18);
Dep_HK_Northeast=rep(0,18);
Int_HK_Northeast=rep(0,18);
Dep_HK_EasternCoastalAreas=rep(0,18);
Int_HK_EasternCoastalAreas=rep(0,18);
Dep_HK_MiddleReachesoftheYangzseRiver=rep(0,18);
Int_HK_MiddleReachesoftheYangzseRiver=rep(0,18);
Dep_HK_SouthernCoastalAreas=rep(0,18);
Int_HK_SouthernCoastalAreas=rep(0,18);
Dep_HK_SoutheastZone=rep(0,18);
Int_HK_SoutheastZone=rep(0,18);
Dep_HK_NorthwestZone=rep(0,18);
Int_HK_NorthwestZone=rep(0,18);


for(t in 1:18)
{
dataList_MAIN_TW_NorthernCoastalAreas=dataList_MAIN_TW_HK[[t]];
dataList_MAIN_HK_NorthernCoastalAreas=dataList_MAIN_TW_HK[[t]];

DD_TW_NorthernCoastalAreas=dataList_MAIN_TW_NorthernCoastalAreas[NorthernCoastalAreas_TW,];
DD_HK_NorthernCoastalAreas=dataList_MAIN_HK_NorthernCoastalAreas[NorthernCoastalAreas_HK,];
DD_TW_MiddleReachesoftheYellowRiver=dataList_MAIN_TW_NorthernCoastalAreas[MiddleReachesoftheYellowRiver_TW,];
DD_HK_MiddleReachesoftheYellowRiver=dataList_MAIN_HK_NorthernCoastalAreas[MiddleReachesoftheYellowRiver_HK,];
DD_TW_Northeast=dataList_MAIN_TW_NorthernCoastalAreas[Northeast_TW,];
DD_HK_Northeast=dataList_MAIN_HK_NorthernCoastalAreas[Northeast_HK,];
DD_TW_EasternCoastalAreas=dataList_MAIN_TW_NorthernCoastalAreas[EasternCoastalAreas_TW,];
DD_HK_EasternCoastalAreas=dataList_MAIN_HK_NorthernCoastalAreas[EasternCoastalAreas_HK,];
DD_TW_MiddleReachesoftheYangzseRiver=dataList_MAIN_TW_NorthernCoastalAreas[MiddleReachesoftheYangzseRiver_TW,];
DD_HK_MiddleReachesoftheYangzseRiver=dataList_MAIN_HK_NorthernCoastalAreas[MiddleReachesoftheYangzseRiver_HK,];
DD_TW_SouthernCoastalAreas=dataList_MAIN_TW_NorthernCoastalAreas[SouthernCoastalAreas_TW,];
DD_HK_SouthernCoastalAreas=dataList_MAIN_HK_NorthernCoastalAreas[SouthernCoastalAreas_HK,];
DD_TW_SoutheastZone=dataList_MAIN_TW_NorthernCoastalAreas[SoutheastZone_TW,];
DD_HK_SoutheastZone=dataList_MAIN_HK_NorthernCoastalAreas[SoutheastZone_HK,];
DD_TW_NorthwestZone=dataList_MAIN_TW_NorthernCoastalAreas[NorthwestZone_TW,];
DD_HK_NorthwestZone=dataList_MAIN_HK_NorthernCoastalAreas[NorthwestZone_HK,];

Dep_TW_NorthernCoastalAreas[t]=DepInt(DD_TW_NorthernCoastalAreas,dim(DD_TW_NorthernCoastalAreas)[1])[1]; 
Int_TW_NorthernCoastalAreas[t]=DepInt(DD_TW_NorthernCoastalAreas,dim(DD_TW_NorthernCoastalAreas)[1])[2]; 
Dep_TW_MiddleReachesoftheYellowRiver[t]=DepInt(DD_TW_MiddleReachesoftheYellowRiver,dim(DD_TW_MiddleReachesoftheYellowRiver)[1])[1]; 
Int_TW_MiddleReachesoftheYellowRiver[t]=DepInt(DD_TW_MiddleReachesoftheYellowRiver,dim(DD_TW_MiddleReachesoftheYellowRiver)[1])[2]; 
Dep_TW_Northeast[t]=DepInt(DD_TW_Northeast,dim(DD_TW_Northeast)[1])[1]; 
Int_TW_Northeast[t]=DepInt(DD_TW_Northeast,dim(DD_TW_Northeast)[1])[2]; 
Dep_TW_EasternCoastalAreas[t]=DepInt(DD_TW_EasternCoastalAreas,dim(DD_TW_EasternCoastalAreas)[1])[1]; 
Int_TW_EasternCoastalAreas[t]=DepInt(DD_TW_EasternCoastalAreas,dim(DD_TW_EasternCoastalAreas)[1])[2]; 
Dep_TW_MiddleReachesoftheYangzseRiver[t]=DepInt(DD_TW_MiddleReachesoftheYangzseRiver,dim(DD_TW_MiddleReachesoftheYangzseRiver)[1])[1]; 
Int_TW_MiddleReachesoftheYangzseRiver[t]=DepInt(DD_TW_MiddleReachesoftheYangzseRiver,dim(DD_TW_MiddleReachesoftheYangzseRiver)[1])[2]; 
Dep_TW_SouthernCoastalAreas[t]=DepInt(DD_TW_SouthernCoastalAreas,dim(DD_TW_SouthernCoastalAreas)[1])[1]; 
Int_TW_SouthernCoastalAreas[t]=DepInt(DD_TW_SouthernCoastalAreas,dim(DD_TW_SouthernCoastalAreas)[1])[2]; 
Dep_TW_SoutheastZone[t]=DepInt(DD_TW_SoutheastZone,dim(DD_TW_SoutheastZone)[1])[1]; 
Int_TW_SoutheastZone[t]=DepInt(DD_TW_SoutheastZone,dim(DD_TW_SoutheastZone)[1])[2]; 
Dep_TW_NorthwestZone[t]=DepInt(DD_TW_SoutheastZone,dim(DD_TW_SoutheastZone)[1])[1]; 
Int_TW_NorthwestZone[t]=DepInt(DD_TW_SoutheastZone,dim(DD_TW_SoutheastZone)[1])[2]; 


Dep_HK_NorthernCoastalAreas[t]=DepInt(DD_HK_NorthernCoastalAreas,dim(DD_HK_NorthernCoastalAreas)[1])[1]; 
Int_HK_NorthernCoastalAreas[t]=DepInt(DD_HK_NorthernCoastalAreas,dim(DD_HK_NorthernCoastalAreas)[1])[2]; 
Dep_HK_MiddleReachesoftheYellowRiver[t]=DepInt(DD_HK_MiddleReachesoftheYellowRiver,dim(DD_HK_MiddleReachesoftheYellowRiver)[1])[1]; 
Int_HK_MiddleReachesoftheYellowRiver[t]=DepInt(DD_HK_MiddleReachesoftheYellowRiver,dim(DD_HK_MiddleReachesoftheYellowRiver)[1])[2]; 
Dep_HK_Northeast[t]=DepInt(DD_HK_Northeast,dim(DD_HK_Northeast)[1])[1]; 
Int_HK_Northeast[t]=DepInt(DD_HK_Northeast,dim(DD_HK_Northeast)[1])[2]; 
Dep_HK_EasternCoastalAreas[t]=DepInt(DD_HK_EasternCoastalAreas,dim(DD_HK_EasternCoastalAreas)[1])[1]; 
Int_HK_EasternCoastalAreas[t]=DepInt(DD_HK_EasternCoastalAreas,dim(DD_HK_EasternCoastalAreas)[1])[2]; 
Dep_HK_MiddleReachesoftheYangzseRiver[t]=DepInt(DD_HK_MiddleReachesoftheYangzseRiver,dim(DD_HK_MiddleReachesoftheYangzseRiver)[1])[1]; 
Int_HK_MiddleReachesoftheYangzseRiver[t]=DepInt(DD_HK_MiddleReachesoftheYangzseRiver,dim(DD_HK_MiddleReachesoftheYangzseRiver)[1])[2]; 
Dep_HK_SouthernCoastalAreas[t]=DepInt(DD_HK_SouthernCoastalAreas,dim(DD_HK_SouthernCoastalAreas)[1])[1]; 
Int_HK_SouthernCoastalAreas[t]=DepInt(DD_HK_SouthernCoastalAreas,dim(DD_HK_SouthernCoastalAreas)[1])[2]; 
Dep_HK_SoutheastZone[t]=DepInt(DD_HK_SoutheastZone,dim(DD_HK_SoutheastZone)[1])[1]; 
Int_HK_SoutheastZone[t]=DepInt(DD_HK_SoutheastZone,dim(DD_HK_SoutheastZone)[1])[2]; 
Dep_HK_NorthwestZone[t]=DepInt(DD_HK_SoutheastZone,dim(DD_HK_SoutheastZone)[1])[1]; 
Int_HK_NorthwestZone[t]=DepInt(DD_HK_SoutheastZone,dim(DD_HK_SoutheastZone)[1])[2]; 
}


par(mfrow=c(2,2))
scale=2003:2020;
plot(scale,Dep_TW_NorthernCoastalAreas,    lwd=2,type="b",xlab="Year",ylab="Dependency",col="red"  ,main="TW's dependency on EER",ylim=c(-1,1));
lines(scale,Dep_TW_MiddleReachesoftheYellowRiver,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="blue");
lines(scale,Dep_TW_Northeast,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="green");
lines(scale,Dep_TW_EasternCoastalAreas,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="yellow");
lines(scale,Dep_TW_MiddleReachesoftheYangzseRiver,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="violet");
lines(scale,Dep_TW_SouthernCoastalAreas,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="black");
lines(scale,Dep_TW_SoutheastZone,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="violet");
lines(scale,Dep_TW_NorthwestZone,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="black");


plot(scale,Dep_HK_NorthernCoastalAreas,    lwd=2,type="b",xlab="Year",ylab="Dependency",col="red"  ,main="HK's dependency on EER",ylim=c(-1,1));
lines(scale,Dep_HK_MiddleReachesoftheYellowRiver,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="blue");
lines(scale,Dep_HK_Northeast,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="green");
lines(scale,Dep_HK_EasternCoastalAreas,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="yellow");
lines(scale,Dep_HK_MiddleReachesoftheYangzseRiver,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="violet");
lines(scale,Dep_HK_SouthernCoastalAreas,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="black");
lines(scale,Dep_HK_SoutheastZone,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="violet");
lines(scale,Dep_HK_NorthwestZone,   lwd=2,type="b",xlab="Year",ylab="Dependency",col="black");

plot(scale,Int_TW_NorthernCoastalAreas,    lwd=2,type="b",xlab="Year",ylab="Integration",col="red" ,main="TW's integration into EER",ylim=c(-1,1));
lines(scale,Int_TW_MiddleReachesoftheYellowRiver,   lwd=2,type="b",xlab="Year",ylab="Integration",col="blue");
lines(scale,Int_TW_Northeast,   lwd=2,type="b",xlab="Year",ylab="Integration",col="green");
lines(scale,Int_TW_EasternCoastalAreas,   lwd=2,type="b",xlab="Year",ylab="Integration",col="yellow");
lines(scale,Int_TW_MiddleReachesoftheYangzseRiver,   lwd=2,type="b",xlab="Year",ylab="Integration",col="violet");
lines(scale,Int_TW_SouthernCoastalAreas,   lwd=2,type="b",xlab="Year",ylab="Integration",col="black");
lines(scale,Int_TW_SoutheastZone,   lwd=2,type="b",xlab="Year",ylab="Integration",col="violet");
lines(scale,Int_TW_NorthwestZone,   lwd=2,type="b",xlab="Year",ylab="Integration",col="black");

plot(scale,Int_HK_NorthernCoastalAreas,    lwd=2,type="b",xlab="Year",ylab="Integration",col="red" ,main="HK dep. Northern Coastal Areas",ylim=c(-1,1));
lines(scale,Int_HK_MiddleReachesoftheYellowRiver,   lwd=2,type="b",xlab="Year",ylab="Integration",col="blue" ,main="HK's integration into EER");
lines(scale,Int_HK_Northeast,   lwd=2,type="b",xlab="Year",ylab="Integration",col="green");
lines(scale,Int_HK_EasternCoastalAreas,   lwd=2,type="b",xlab="Year",ylab="Integration",col="yellow");
lines(scale,Int_HK_MiddleReachesoftheYangzseRiver,   lwd=2,type="b",xlab="Year",ylab="Integration",col="violet");
lines(scale,Int_HK_SouthernCoastalAreas,   lwd=2,type="b",xlab="Year",ylab="Integration",col="black");
lines(scale,Int_HK_SoutheastZone,   lwd=2,type="b",xlab="Year",ylab="Integration",col="violet");
lines(scale,Int_HK_NorthwestZone,   lwd=2,type="b",xlab="Year",ylab="Integration",col="black");


















EER=data1_MAIN_TW_HK[1:31,6];
NorthernCoastalAreas=which(EER=="Northern Coastal Areas");
MiddleReachesoftheYellowRiver=which(EER=="Middle Reaches of the Yellow River");
Northeast=which(EER=="Northeast ");
EasternCoastalAreas=which(EER=="Eastern Coastal Areas");
SouthernCoastalAreas=which(EER=="Southern Coastal Areas");
MiddleReachesoftheYangzseRiver=which(EER=="Middle Reaches of the Yangzse River");
SoutheastZone=which(EER=="Southeast Zone");
NorthwestZone=which(EER=="Northwest Zone");


##---categorical dependency and integration

TER=data1_MAIN_TW_HK[1:31,3];
Eastern=which(TER=="Eastern Regions");




%%par(mfrow=c(4,1));
scale1=2004:2020;
plot(scale1,vec_MAIN_TW_HK,    lwd=2,type="b",xlab="Year",ylab="Feature distance",col="red"  ,main="Consecutive feature shift",ylim=c(0,2.5));
lines(scale1,vec_MAIN_TW,      lwd=2,type="b",xlab="Year",ylab="Feature distance",col="blue" ,main="Consecutive feature shift");
lines(scale1,vec_MAIN_HK,      lwd=2,type="b",xlab="Year",ylab="Feature distance",col="green",main="Consecutive feature shift");
lines(scale1,vec_MAIN,         lwd=2,type="b",xlab="Year",ylab="Feature distance",col="black",main="Consecutive feature shift");
legend("topleft",legend=c("Mainland with TW and HK","Mainland with TW","Mainland with HK","Mainland"), 
col=c("red", "blue","green","black"),lwd=2,lty=1:2, cex=0.8);



len=length(PcaList_MAIN_TW_HK);
len_MAIN_TW=length(PcaList_MAIN_TW);
len_MAIN_HK=length(PcaList_MAIN_HK);
len_MAIN=length(PcaList_MAIN);

dis12=dis13=dis14=dis23=dis24=dis34=rep(0,len);  ##1:MAIN_TW_HK, 2:MAIN_TW, 3:MAIN_HK, 4:MAIN

for(i in 1:len_MAIN_TW_HK)
{
PCA_1=PcaList_MAIN_TW_HK[[i]];
PCA_2=PcaList_MAIN_TW[[i]];
PCA_3=PcaList_MAIN_HK[[i]];
PCA_4=PcaList_MAIN[[i]];

dist12=FeatureDis(PCA_1,PCA_2);
dist13=FeatureDis(PCA_1,PCA_3);
dist14=FeatureDis(PCA_1,PCA_4);
dist23=FeatureDis(PCA_2,PCA_3);
dist24=FeatureDis(PCA_2,PCA_4);
dist34=FeatureDis(PCA_3,PCA_4);
 
dis12[i]=dist12;
dis13[i]=dist13;
dis14[i]=dist14;
dis23[i]=dist23;
dis24[i]=dist24;
dis34[i]=dist34;
}

scale2=2003:2020;
par(mfrow=c(2,3));
plot(scale2,dis12, main="Annaul distance between F1,F2",type="b",xlab="Year",ylab="Distance between features");
plot(scale2,dis13,  main="Annaul distance between F1,F3",type="b",xlab="Year",ylab="Distance between features");
plot(scale2,dis14,  main="Annaul distance between F1,F4",type="b",xlab="Year",ylab="Distance between features");
plot(scale2,dis23,  main="Annaul distance between F2,F3",type="b",xlab="Year",ylab="Distance between features");
plot(scale2,dis24,  main="Annaul distance between F2,F4",type="b",xlab="Year",ylab="Distance between features");
plot(scale2,dis34,  main="Annaul distance between F3,F4",type="b",xlab="Year",ylab="Distance between features");



timeSeries1=t(apply(data111,1,rev));
timeSeries2=t(apply(data222,1,rev));
timeSeries3=t(apply(data333,1,rev));

TER=as.list(numeric(3));dim(TER)=c(3,1);
Eastern=data1[,3]=="Eastern Regions";
IND1=which(Eastern==TRUE);
Central=data1[,3]=="Central Regions";
IND2=which(Central==TRUE);
Western=data1[,3]=="Western Regions";
IND3=which(Western==TRUE);
TER[[1,1]]=IND1;
TER[[2,1]]=IND2;
TER[[3,1]]=IND3;

n1=length(IND1);
n2=length(IND2);
n3=length(IND3);

TimeSeries_Eastern=as.list(numeric(20));dim(TimeSeries_Eastern)=c(1,20);
TimeSeries_Central=as.list(numeric(20));dim(TimeSeries_Central)=c(1,20);
TimeSeries_Western=as.list(numeric(20));dim(TimeSeries_Western)=c(1,20);


FeatureList=as.list(numeric(31*1)); dim(FeatureList)=c(31,1);
for(i in 1:31)
{
mat=matrix(nrow=3,ncol=20);
mat[1,]=timeSeries1[i,];
mat[2,]=timeSeries2[i,];
mat[3,]=timeSeries3[i,];
FeatureList[[i,1]]=mat;
}

FeatureList1=FeatureList[TER[[1,1]]];
FeatureList2=FeatureList[TER[[2,1]]];
FeatureList3=FeatureList[TER[[3,1]]];

for(i in 1:20)
{
mat_east=matrix(nrow=3,ncol=11);
mat_centre=matrix(nrow=3,ncol=8);
mat_west=matrix(nrow=3,ncol=12);
   for(j in 1:11)
   {
   mat_east[,j]=FeatureList1[[j]][,i]
   }
   for(j in 1:8)
   {
   mat_centre[,j]=FeatureList2[[j]][,i]
   }
   for(j in 1:12)
   {
   mat_west[,j]=FeatureList3[[j]][,i]
   }
TimeSeries_Eastern[[i]]=mat_east;
TimeSeries_Central[[i]]=mat_centre;
TimeSeries_Western[[i]]=mat_west;
}

scatterplot3d(t(TimeSeries_Eastern[[1]]),pch=16,color="red",box=FALSE);
scatterplot3d(t(TimeSeries_Central[[1]]),pch=16,color="blue",box=FALSE);
scatterplot3d(t(TimeSeries_Western[[1]]),pch=16,color="black",box=FALSE);

scatterplot3d(t(FeatureList[[1,1]]),ty="o",xlab="Primary Industry",ylab="Secondary Industry",zlab="Tertiary Industry");
n=length(TER[[1,1]]);
layout(matrix(1:3,ncol=3));
scatterplot3d(t(TimeSeries_Eastern[[1]]),pch=16,color="red",box=FALSE);
scatterplot3d(t(TimeSeries_Central[[1]]),pch=16,color="blue",box=FALSE);
scatterplot3d(t(TimeSeries_Western[[1]]),pch=16,color="black",box=FALSE);




for(i in TER[[1,1]])
{
plot3d(t(FeatureList[[i,1]]),ty="o",xlab="Primary Industry",ylab="Secondary Industry",zlab="Tertiary Industry");
}







###printing for primary/secondary/tertiary-industry-against-three-economic-zone GDP 
par(mfrow=c(3,3));
Eastern=data1[,3]=="Eastern Regions";
IND1=which(Eastern==TRUE);
data1_1=data11[IND1,]; 
data1_11=data1_1[,2:21];
data1_111=do.call(cbind,data1_11);
data1_1111=apply(data1_111,1,rev);
data1_11111=t(data1_1111);
plot(2002:2021,  data1_11111[1,],ty="l",main="Primary: Eastern",ylim=c(50,5000),xlab="Year",ylab="GDP in RMB");
points(2002:2021, data1_11111[2,],ty="l");
points(2002:2021, data1_11111[3,],ty="l");
points(2002:2021,data1_11111[4,],ty="l");
points(2002:2021,data1_11111[5,],ty="l");
points(2002:2021,data1_11111[6,],ty="l");
points(2002:2021,data1_11111[7,],ty="l");
points(2002:2021,data1_11111[8,],ty="l");
points(2002:2021,data1_11111[9,],ty="l");


Central=data1[,3]=="Central Regions";
IND2=which(Central==TRUE);
data1_2=data11[IND2,]; 
data1_22=data1_2[,2:21];
data1_222=do.call(cbind,data1_22);
data1_2222=apply(data1_222,1,rev);
data1_2222=t(data1_2222);
plot(2002:2021,data1_2222[1,],ty="l",ylim=c(50,5000),main="Primary: Central",xlab="Year",ylab="GDP in RMB");
points(data1_2222[2,],ty="l");
points(2002:2021,data1_2222[3,],ty="l");
points(2002:2021,data1_2222[4,],ty="l");
points(2002:2021,data1_2222[5,],ty="l");
points(2002:2021,data1_2222[6,],ty="l");
points(2002:2021,data1_2222[7,],ty="l");
points(2002:2021,data1_2222[8,],ty="l");



Western=data1[,3]=="Western Regions";
IND3=which(Western==TRUE);
data1_3=data11[IND3,]; 
data1_33=data1_3[,2:21];
data1_333=do.call(cbind,data1_33);
data1_3333=apply(data1_333,1,rev);
data1_3333=t(data1_3333);
plot(2002:2021,data1_3333[1,],ty="l",ylim=c(50,5000),main="Primary: Western",xlab="Year",ylab="GDP in RMB");
points(2002:2021,data1_3333[2,],ty="l");
points(2002:2021,data1_3333[3,],ty="l");
points(2002:2021,data1_3333[4,],ty="l");
points(2002:2021,data1_3333[5,],ty="l");
points(2002:2021,data1_3333[6,],ty="l");
points(2002:2021,data1_3333[7,],ty="l");
points(2002:2021,data1_3333[8,],ty="l");



Eastern=data2[,3]=="Eastern Regions";
IND1=which(Eastern==TRUE);
data2_1=data22[IND1,]; 
data2_11=data2_1[,2:21];
data2_111=do.call(cbind,data2_11);
data2_1111=apply(data2_111,1,rev);
data2_11111=t(data2_1111);
plot(2002:2021,data2_11111[1,],ty="l",main="Secondary: Eastern",ylim=c(50,50000),xlab="Year",ylab="GDP in RMB");
points(2002:2021,data2_11111[2,],ty="l");
points(2002:2021,data2_11111[3,],ty="l");
points(2002:2021,data2_11111[4,],ty="l");
points(2002:2021,data2_11111[5,],ty="l");
points(2002:2021,data2_11111[6,],ty="l");
points(2002:2021,data2_11111[7,],ty="l");
points(2002:2021,data2_11111[8,],ty="l");
points(2002:2021,data2_11111[9,],ty="l");



Central=data2[,3]=="Central Regions";
IND2=which(Central==TRUE);
data2_2=data22[IND2,]; 
data2_22=data2_2[,2:21];
data2_222=do.call(cbind,data2_22);
data2_2222=apply(data2_222,1,rev);
data2_2222=t(data2_2222);
plot(2002:2021,data2_2222[1,],ty="l",main="Secondary: Central",ylim=c(50,50000),xlab="Year",ylab="GDP in RMB");
points(2002:2021,data2_2222[2,],ty="l");
points(2002:2021,data2_2222[3,],ty="l");
points(2002:2021,data2_2222[4,],ty="l");
points(2002:2021,data2_2222[5,],ty="l");
points(2002:2021,data2_2222[6,],ty="l");
points(2002:2021,data2_2222[7,],ty="l");
points(2002:2021,data2_2222[8,],ty="l");



Western=data2[,3]=="Western Regions";
IND3=which(Western==TRUE);
data2_3=data22[IND3,]; 
data2_33=data2_3[,2:21];
data2_333=do.call(cbind,data2_33);
data2_3333=apply(data2_333,1,rev);
data2_3333=t(data2_3333);
plot(2002:2021,data2_3333[1,],ty="l",main="Secondary: Western",ylim=c(50,50000),xlab="Year",ylab="GDP in RMB");
points(2002:2021,data2_2222[2,],ty="l");
points(2002:2021,data2_2222[3,],ty="l");
points(2002:2021,data2_2222[4,],ty="l");
points(2002:2021,data2_2222[5,],ty="l");
points(2002:2021,data2_2222[6,],ty="l");
points(2002:2021,data2_2222[7,],ty="l");
points(2002:2021,data2_2222[8,],ty="l");
points(2002:2021,data2_3333[2,],ty="l");
points(2002:2021,data2_3333[3,],ty="l");
points(2002:2021,data2_3333[4,],ty="l");
points(2002:2021,data2_3333[5,],ty="l");
points(2002:2021,data2_3333[6,],ty="l");
points(2002:2021,data2_3333[7,],ty="l");
points(2002:2021,data2_3333[8,],ty="l");




Eastern=data3[,3]=="Eastern Regions";
IND1=which(Eastern==TRUE);
data3_1=data33[IND1,]; 
data3_11=data3_1[,2:21];
data3_111=do.call(cbind,data3_11);
data3_1111=apply(data3_111,1,rev);
data3_11111=t(data3_1111);
plot(2002:2021,data3_11111[1,],ty="l",main="Tertiary: Eastern",ylim=c(50,50000),xlab="Year",ylab="GDP in RMB");
points(2002:2021,data3_11111[2,],ty="l");
points(2002:2021,data3_11111[3,],ty="l");
points(2002:2021,data3_11111[4,],ty="l");
points(2002:2021,data3_11111[5,],ty="l");
points(2002:2021,data3_11111[6,],ty="l");
points(2002:2021,data3_11111[7,],ty="l");
points(2002:2021,data3_11111[8,],ty="l");
points(2002:2021,data3_11111[9,],ty="l");



Central=data3[,3]=="Central Regions";
IND2=which(Central==TRUE);
data3_2=data33[IND2,]; 
data3_22=data3_2[,2:21];
data3_222=do.call(cbind,data3_22);
data3_2222=apply(data3_222,1,rev);
data3_2222=t(data3_2222);
plot(2002:2021,data3_2222[1,],ty="l",main="Tertiary: Central",ylim=c(50,50000),xlab="Year",ylab="GDP in RMB");
points(2002:2021,data3_2222[2,],ty="l");
points(2002:2021,data3_2222[3,],ty="l");
points(2002:2021,data3_2222[4,],ty="l");
points(2002:2021,data3_2222[5,],ty="l");
points(2002:2021,data3_2222[6,],ty="l");
points(2002:2021,data3_2222[7,],ty="l");
points(2002:2021,data3_2222[8,],ty="l");



Western=data3[,3]=="Western Regions";
IND3=which(Western==TRUE);
data3_3=data33[IND3,]; 
data3_33=data3_3[,2:21];
data3_333=do.call(cbind,data3_33);
data3_3333=apply(data3_333,1,rev);
data3_3333=t(data3_3333);
plot(2002:2021,data3_3333[1,],ty="l",main="Tertiary: Western",ylim=c(50,50000),xlab="Year",ylab="GDP in RMB");
points(2002:2021,data3_3333[2,],ty="l");
points(2002:2021,data3_3333[3,],ty="l");
points(2002:2021,data3_3333[4,],ty="l");
points(2002:2021,data3_3333[5,],ty="l");
points(2002:2021,data3_3333[6,],ty="l");
points(2002:2021,data3_3333[7,],ty="l");
points(2002:2021,data3_3333[8,],ty="l");






































##trial
Dmat=matrix(sample(-129:43,50),10,5);
clVec=matrix(sample(1:5,5),5,1);
MahaMetric(clVec,DMat,ty="l");
 
##trial
B1=matrix(sample(-100:100,25),5,5);
B2=matrix(sample(-100:100,25),5,5);
B3=matrix(sample(-100:100,25),5,5);

p_colMat=matrix(sample(-12:23,5),5,1);
q_colMat=matrix(sample(-12:23,5),5,1);
r_colMat=matrix(sample(-12:3,5),5,1);

delta_pq=Delta1(p_colMat,B1,q_colMat,B2);
delta_qr=Delta1(q_colMat,B2,r_colMat,B3);
delta_pr=Delta1(p_colMat,B1,r_colMat,B3);
delta_pq;
delta_qr;
delta_pr;
##trial
Dmat=matrix(sample(-129:43,50),10,5);
Emat=matrix(sample(-129:43,30),6,5);
Fmat=matrix(sample(-129:43,80),16,5);
disMat_DE=DisMat(Dmat,Emat,B1,B2); 
disMat_EF=DisMat(Emat,Fmat,B2,B3); 
disMat_DF=DisMat(Dmat,Fmat,B1,B3); 
##trial
hausdorff_DE=Hausdorff(disMat_DE);
hausdorff_EF=Hausdorff(disMat_EF);
hausdorff_DF=Hausdorff(disMat_DF);
hausdorff_DE;
hausdorff_EF;
hausdorff_DF;
##trial
B1=matrix(sample(-100:100,25),5,5);
B2=matrix(sample(-100:100,25),5,5);
Dmat=matrix(sample(-129:43,50),10,5);
Emat=matrix(sample(-129:43,30),6,5);
delta2=Delta2(Dmat,Emat,B1,B2);
delta2;
 
##trial
alpha=0.4; beta=0.7;
delta3_DE=Delta3(p_colMat,Dmat,q_colMat,Emat,alpha,beta) 
delta3_EF=Delta3(p_colMat,Emat,q_colMat,Fmat,alpha,beta) 
delta3_DF=Delta3(p_colMat,Dmat,q_colMat,Fmat,alpha,beta) 
delta3_DE;
delta3_EF;
delta3_DF;


###



M=matrix(sample(1:3,10*4,replace=TRUE),10,4);
DM=list(M);
DM=data.frame(DM);

color=c("red","green","blue");
colors = color[as.numeric(DM$X4)];
s3d=scatterplot3d(DM[,1:3],pch = 16,color=colors);
legend("bottom", legend = unique(DM$X4), pch = 16,
col=color,inset = -0.25, xpd = TRUE, horiz = TRUE);



s3d =scatterplot3d(iris[,1:3], pch = 16, color)
legend("right", legend = levels(iris$Species),
  color, pch = 16, inset = 0.1)