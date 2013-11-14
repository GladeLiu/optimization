revisenewton<-function(fun,gfun,Hess,x0,maxk=150,epsilon=1e-5){
#Revise Newton Method
  n<-length(x0); rho<-0.55;
  sigma<-0.4; tau<-0.0;k<-0;
  while(k<maxk){
    gk<-gfun(x0);
    muk<-norm(gk)^(1+tau);
    Gk<-Hess(x0);
    Ak<-Gk+muk*diag(n);
    dk<-solve(Gk,-gk); 
    if(norm(gk)<epsilon){
       break; 
    }
    m<-0;mk<-0;
    while(m<20){   #Armijo search 
        if(fun(x0+rho^m*dk)<fun(x0)+sigma*rho^m*t(gk)%*%dk){
            mk<-m;break;
        }
        m<-m+1;
    }
    x0<-x0+rho^mk*dk;
    k<-k+1;
  }
  x<-x0;
  val<-fun(x); 
  list(point=x, value=val);
}

#A simple example
#objective function
fun<-function(x){
  f<-100*(x[1]^2-x[2])^2+(x[1]-1)^2;
}
#gradient
gfun<-function(x){
  g<-c(400*x[1]*(x[1]^2-x[2])+2*(x[1]-1), -200*(x[1]^2-x[2]));
  g<-t(t(g));
}
#Hesse Matrix
Hess<-function(x){
  n<-length(x);
  He<-matrix(0,n,n);
  He<-matrix(c(1200*x[1]^2-400*x[2]+2,-400*x[1],-400*x[1],200),2,byrow=T);
}
#initial value
x0<-c(1.4,0.4);
revisenewton(fun,gfun,Hess,x0);
