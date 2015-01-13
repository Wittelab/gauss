/* This procedure generates a matrix y of binomial
   variates with probabilities p and totals n. */
proc rndbin(p,n);
local nr,nc,j,k,y;
nr = maxc(rows(p)|rows(n));
nc = maxc(cols(p)|cols(n));
p = ones(nr,nc).*p; n = ones(nr,nc).*n;
y = zeros(nr,nc);
j=0;
do until j ge nr; j = j + 1;
k = 0;
     do until k ge nc; k = k + 1;
          if n[j,k]; y[j,k] = sumc(rndu(n[j,k],1) .lt p[j,k]); endif;
     endo;
endo;
retp(y);
endp;