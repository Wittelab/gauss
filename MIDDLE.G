/* This proc finds the middle of the three values in each row of x~y~z,
   where x,y, & z are vectors. Runs in under 7% the time of median(x). */
proc middle(x,y,z);
local v; v = x'|y'|z';
retp(sumc(v) - maxc(v) - minc(v));
endp;