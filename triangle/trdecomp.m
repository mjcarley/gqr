dat = load("triangle.dat") ;

xt = dat(1:3,1:2) ;
x0 = dat(4, 1:2) ;

idx = [1 2 3 1] ;

r  = dat(5:end,1) ;
th = dat(5:end,2) ;
r0 = dat(5:end,3) ;
th0 = dat(5:end,4) ;

z0 = x0*[1; j] ;

zt = [] ;
for i=1:3
  zt = [zt [0; r(i)*exp(j*th(i)); r(i)*r0(i)*exp(j*(th0(i)+th(i)))]] ;
  ##zt = [zt [0; r(i)*exp(j*th(i)); r(i)*r0(i)*exp(j*(th0(i)))]] ;
endfor

zt += z0 ;
