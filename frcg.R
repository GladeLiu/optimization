frcg<-function(fun,gfun,x0,epsilon=1e-4,maxk=5000){
#FRCG method
  rho<-0.6;sigma<-0.4;
  k<-0;n<-length(x0);
  while(k<maxk){
    g=gfun(x0);
    itern=k-(n+1)*floor(k/(n+1));
    itern=itern+1;
    if(itern==1){  
        d=-g;  
    }else{
        beta=(t(g)%*%g)/(t(g0)%*%g0);
        d=-g+as.vector(beta)*d0;  gd=t(g)%*%d;
        if(gd>=0.0){
            d=-g;  
        }
    }
    if(norm(g)<epsilon){
        break;
    }
    m<-0; mk<-0;
    while(m<20){
        if(fun(x0+rho^m*d)<fun(x0)+sigma*rho^m*t(g)%*%d){
            mk=m; break;
        }
        m<-m+1;
    }
    x0<-x0+rho^mk*d;
    val<-fun(x0);
    g0<-g;d0<-d; 
    k<-k+1;
  }
  x<-x0;
  val<-fun(x0); 
  list(point=x,value=val);
}

#example
#objective function
fun<-function(x){
 f<-100*(x[1]^2-x[2])^2+(x[1]-1)^2
}
#gradient
gfun<-function(x){
 g<-c(400*x[1]*(x[1]^2-x[2])+2*(x[1]-1), -200*(x[1]^2-x[2]));
 g<-t(t(g));
}
#the initial value
x0<-c(1.2,0.6);
frcg(fun,gfun,x0);
