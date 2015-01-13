/* prints a character matrix in row-major 8-column matrix form. */
proc (0) = printnam(n);
local i; 
format /ldn 10,0; n = vecr(n); i=0;
do until i ge rows(n); i=i+1; 
   $n[i];; if not i%8 or i eq rows n; print; endif; 
endo; 
format /rds 9,3;
endp;
