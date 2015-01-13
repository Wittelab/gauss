/* This proc computes values of the logistic-beta density f(x;a,b)
   where a,b are the beta parameters.
   Inputs: x; and a,b = beta parameters
*/
proc lgstbeta(x,a,b); local c;
if sumall((a~b) .le 1);
   c = ln(gamma(a+b)./(gamma(a).*gamma(b)));
   else; c = lnfact(a+b-1)-lnfact(a-1)-lnfact(b-1);
endif;
retp(exp(c + a.*x - (a+b).*ln(1+exp(x))));
endp;
