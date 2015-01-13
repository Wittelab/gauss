/* density of generalized ed log chi-squared variate x with df degrees of freedom, 
   location m, scale s.
   equals exp(-r)*(r^a)/a!s where r = exp((x-m)/s)+ln(a)) and a = df/2 */
proc logchi2(x,df,m,s);
local a,g; a = df/2;
g = (x-m)./s+ln(a);
retp(exp(a.*g-exp(g)-lngamma(a))./s);
endp;