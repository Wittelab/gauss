/* This proc returns the generalized logit transform of the matrix p
   of nonzero proportions. Each row of p must sum to less than 1.
   (useful for multinomial logistic modelling) */
/* Test expitg and logitg: */
a = seqa(1,1,5)/5; p = expitg(a,1);
format /rdn 10,3; a~expit(a)~p~logit(p)~logitg(p);
p = expitg(a~(5*a),ones(2,1)); a~p~logitg(p);
end;
proc logitg(p);
local s;
s = sumr(p);
if sumc(s .ge 1) or sumall(p .le 0);
   "ERROR: Arguments out of range in logitp.g"; end; endif;
retp(ln(p./(1-s)));
endp;

