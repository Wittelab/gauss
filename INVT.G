/* Test program:
let p = .025 .5 .95 .975 .99 .995 .995 .975;
let df = 8    4   3   3   3    4     3   2;
t = invt(p,df); p~df~t~(1-cdftc(t,df)); end; */
/* This proc inverts a t distribution by brute search.
   Accurate to at least two digits past decimal point
   for df ge 3 and .01 < p < .99.
*/
proc (1) = invt(p,df);
local signt,z,t;
signt = p .ge .5;
p = signt.*p + (1-signt).*(1-p);
signt = 2*signt - 1;
if sumc(p .gt .99);
   "WARNING: INVT.G can be inaccurate for p > .99 or p < .01"; endif;
if sumc(df .lt 3);
   "WARNING: INVT.G is inaccurate for df < 3"; endif;
z = invnorm(p) + seqa(0,.0049,500)';
t = 1-cdftc(z,df.*ones(1,500));
t = sumr(t .lt p);
retp(signt.*(z[.,1] + .0049*t));
endp;
