function dI=tfunc(t)

global x y n

R = sqrt((t-x).^2 + y.^2) ;
dI = t.^n./R.^3 ;

