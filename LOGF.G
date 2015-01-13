/* density of generalized log F variate x with dfa,dfb degrees of freedom, 
   location m, scale s.
   equals [ab/(a+b)]C(a+b,a)(p^a)(1-p)^b where p = expit((x-m)/s+ln(a/b)). */
proc logF(x,dfa,dfb,m,s);
local a,b,p; a = dfa/2; b = dfb/2; 
p = expit((x-m)./s+ln(a./b));
retp(exp(a.*ln(p) + b.*ln(1-p) - ln(beta(a,b)))./s);
endp;