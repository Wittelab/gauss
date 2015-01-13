/* This proc computes values of the beta density f(p;a,b) at p
   where a,b are the beta parameters.
*/
proc densbeta(p,a,b);
if sumall(p .lt 0) or sumall(p .gt 1);"INPUT ERROR: p < 0 OR p > 1"; end; endif;
a = a - 1; b = b - 1;
retp(exp(lnfact(a+b+1) - lnfact(a) - lnfact(b) + a.*ln(p) + b.*ln(1-p)));
endp;
