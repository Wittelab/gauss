/* This proc solves the quadratic equation a*x^2 + b*x + c = 0.
   Returns the roots (which may be complex) in a two-column matrix
   (one row for each equation).
*/
proc quadsoln(a,b,c);
local a0,b0,d;
a0 = a .eq 0; b0 = b .eq 0;
if sumc(a0 .and b0);
   "ERROR: a and b cannot both be zero"; end; endif;
d = sqrt(b.^2 - 4*a.*c); d = d~(-d);
retp((((1-a0).*(-b + d) - a0.*c./(b + b0)))./(2*a + a0));
endp;
