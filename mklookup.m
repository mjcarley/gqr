n = 1024 ;

th = linspace(-pi-0.5, pi+0.5, n) ;
f = sin(th) ;

fid = fopen('tmp', 'w') ;

fprintf(fid, 'const gdouble _GQR_SINE_LOOKUP[] = \n  {') ;
fprintf(fid, '%1.20e, %1.20e,\n', f(1:2)) ;
fprintf(fid, '   %1.20e, %1.20e,\n', f(3:end-1)) ;
fprintf(fid, '%1.20e} ;\n', f(end)) ;

fprintf(fid, 
	'#define _GQR_SINE_1_DELTA_THETA %1.20e\n',
	1/diff(th(1:2))) ;
fprintf(fid, 
	'#define _GQR_SINE_DELTA_THETA %1.20e\n',
	diff(th(1:2))) ;
fprintf(fid, 
	'#define _GQR_SINE_THETA_MIN %1.20e\n', th(1)) ;

fprintf(fid, 
	'#define _GQR_SINE_THETA_MAX %1.20e\n', th(end)) ;


fclose(fid) ;


