function H=hermite(n,x)

x = x(:)' ;

H = ones(n+1,length(x)) ;

H(2,:) = 2*x ;

for i=1:n-1
  H(i+2,:) = 2*x.*H(i+1,:) - 2*i*H(i,:) ;
end