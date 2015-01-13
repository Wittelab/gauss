/* computes logistic density for p with logit location m, scale s */ 
proc logistic(p,m,s);
local f;
f = expit((logit(p)-m)./s); 
retp(f.*(1-f)./(s.*p.*(1-p)));
endp; 
