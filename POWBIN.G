/* Test (should yield .95,.95):  powbin(.25,.4,.01,530|357,265|357)'; */
/* This procedure calculates power of a two-group two-sided binomial
   test, based on solving the approximate formulas 3.19-3.20 of Fleiss (1981).
   Inputs: p1, p0 = outcome probabilities in the two groups under alternative;
           alpha = desired alpha-level;
           n1,n0 = sample sizes of each group.
   Outputs: power.
*/
proc powbin(p1,p0,alpha,n1,n0);
local ca,r,da,pm,cb;
ca = cinvnorm(alpha/2);
r = n0./n1;  da = abs(p1 - p0);
n1 = n1 - (1+r)./(r.*da);
pm = (p1 + r.*p0)./(1 + r);
cb = ca.*sqrt((1+r).*pm.*(1-pm)) - sqrt(n1.*r).*da;
cb = cb./sqrt(r.*p1.*(1-p1) + p0.*(1-p0));
retp(cdfnc(cb));
endp;
