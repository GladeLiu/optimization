grad<-function(fun,gfun,x0,epsilon=1e-10,maxk<-5000){ 
  rho<-0.5;sigma<-0.4;k<-0; 
  while(k<maxk){
    g=gfun(x0); #calculating grdient
    d=-g;    
    if(norm(d)<epsilon){
      break;
    }
    m<-0; mk<-0;
    while(m<20){   #Armijo search
        if(fun(x0+rho^m*d)<fun(x0)+sigma*rho^m*t(g)%*%d){
            mk<-m; break;
        }
      m<-m+1;
    }
    x0<-x0+rho^mk*d;
    k<-k+1;
  }
  x<-x0;
  val<-fun(x0); 
  list(point=x, value=val);
}

#example
fun<-function(x){
 f<-100*(x[1]^2-x[2])^2+(x[1]-1)^2
}

gfun<-function(x){
 g<-c(400*x[1]*(x[1]^2-x[2])+2*(x[1]-1), -200*(x[1]^2-x[2]));
 g<-t(t(g));
}
grad(fun,gfun,x0);

