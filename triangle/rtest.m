r0 = 0.9 ; th0 = 75*pi/180 ;

z = linspace(1, r0*exp(j*th0), 16).' ;

th = atan2(imag(z), real(z)) ;

M = r0*sin(th0)./(r0*sin(th0 - th) + sin(th)) ;



