/* This procedure generates a rows(p) by r matrix a of r multinomial variates
   with distributions p[.,k] and totals n[k]. 
   Faster than rndmult2 when n < 55. */
proc rndmultn(r,p,n);
local c,a,k,u;
n = n.*ones(r,1);
p = p.*ones(rows(p),r);
c = cumsumc(p);
a = zeros(rows(p),r);
k=0;
do until k ge r; k=k+1;
      u = rndu(n[k],1);
      a[.,k] = sumc((u .gt (c[.,k] - p[.,k])').*(u .le c[.,k]'));
endo;
retp(a);
endp;

