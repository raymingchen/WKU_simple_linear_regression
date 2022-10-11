#### WKU: Simple Linear Regression (SLR)
rm(list=ls());


##Mother nature/data generator (deterministic)+ random part
x=seq(0,19,by=1);
y_determ=1+1.45*x;   
err_rand=rnorm(length(x),0,2);  ## assume the residual distribution is standard normal dist.
y=y_determ+err_rand;

n=length(x);
sum_x=sum(x);
sum_y=sum(y);
sum_xy=x%*%y;
sum_x_2=sum(x^2);

numeratorMat_alpha=matrix(c(sum_y,sum_xy,sum_x,sum_x_2),2,2);
numeratorMat_bata=matrix(c(n,sum_x,sum_y,sum_xy),2,2);
denominatorMat=matrix(c(n,sum_x,sum_x,sum_x_2),2,2);
alpha_hat=det(numeratorMat_alpha)/det(denominatorMat);
beta_hat=det(numeratorMat_bata)/det(denominatorMat);

y_hat=alpha_hat+beta_hat*x;

par(mfrow=c(3,4));
plot(x,y_determ,main="y_determ=f(x)=1+1.45x");
plot(x,err_rand,xlab="time",main="err~N(0,2)");
plot(x,y,main="y[x_i]=y_determ[x_i]+err[i]");
plot(x,y_hat,main="y_hat=alpha_hat+beta_hat*x");
##-------------
n=length(x);
sum_x=sum(x);
sum_y=sum(y);
sum_xy=x%*%y;
sum_x_2=sum(x^2);

numeratorMat_alpha=matrix(c(sum_y,sum_xy,sum_x,sum_x_2),2,2);
numeratorMat_bata=matrix(c(n,sum_x,sum_y,sum_xy),2,2);
denominatorMat=matrix(c(n,sum_x,sum_x,sum_x_2),2,2);
alpha_hat=det(numeratorMat_alpha)/det(denominatorMat);
beta_hat=det(numeratorMat_bata)/det(denominatorMat);


err_rand=rnorm(n,0,2);  ## assume the residual distribution is standard normal dist.
y=y_determ+err_rand;
plot(x,y_determ,main="y_determ=f(x)=1+1.45x");
plot(x,err_rand,xlab="time",main="err~N(0,2)");
plot(x,y,main="y[x_i]=y_determ[x_i]+err[i]");
plot(x,y_hat,main="y_hat=alpha_hat+beta_hat*x");

##-------------
n=length(x);
sum_x=sum(x);
sum_y=sum(y);
sum_xy=x%*%y;
sum_x_2=sum(x^2);

numeratorMat_alpha=matrix(c(sum_y,sum_xy,sum_x,sum_x_2),2,2);
numeratorMat_bata=matrix(c(n,sum_x,sum_y,sum_xy),2,2);
denominatorMat=matrix(c(n,sum_x,sum_x,sum_x_2),2,2);
alpha_hat=det(numeratorMat_alpha)/det(denominatorMat);
beta_hat=det(numeratorMat_bata)/det(denominatorMat);


err_rand=rnorm(n,0,2);  ## assume the residual distribution is standard normal dist.
y=y_determ+err_rand;
plot(x,y_determ,main="y_determ=f(x)=1+1.45x");
plot(x,err_rand,xlab="time",main="err~N(0,2)");
plot(x,y,main="y[x_i]=y_determ[x_i]+err[i]");
plot(x,y_hat,main="y_hat=alpha_hat+beta_hat*x");
