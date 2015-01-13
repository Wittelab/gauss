proc orcrude(a);
local z,r,se;
let z[1,3] = 0 -1.96 1.96;
r = a[.,1].*a[.,4]./(a[.,2].*a[.,3]);
se = sqrt(sumr(1/a));
retp(r.*exp(z.*se));
endp;