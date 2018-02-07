global x y n

N = 6 ;

x = 0.9 ; y = 0.9 ;

I = zeros(N+1,1) ;

for n=0:N
  I(n+1) = quad("tfunc", -1, 1, [1e-9 1e-9]) ;
end