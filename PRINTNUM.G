/* prints a numeric matrix in row-major 8-column matrix form. */
proc (0) = printnum(n);
local i; 
format /rds 10,3; n = vecr(n); i=0;
do until i ge rows(n); i=i+1; 
   n[i];; if not i%8 or i eq rows n; print; endif; 
endo; 
format 9,3;
endp;