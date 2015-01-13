/* Test (should yield 530,357):  ssbin(.25,.4,.01,.95,.5|1)'; */
/* This procedure calculates sample sizes required for a two-group binomial
   test according to the approximate formulas 3.18-3.19 of Fleiss (1981).
   Inputs: p1, p0 = outcome probabilities in the two groups under alternative;
           alpha = desired alpha-level;
           power = desired power;
           r = desired ratio of sample 1 to sample 0.
   Outputs: n1 = sample size required for group 1
                 (sample size for group 0 is r.*n1.
*/
proc ssbin(p1,p0,alpha,power,r);
local ca,cb,pm,da,s,n1;
ca = cinvnorm(alpha/2); cb = cinvnorm(power);
da = abs(p1 - p0);  pm = (p1 + r*p0)./(1 + r);
s = (ca.*sqrt((1+r).*pm.*(1-pm)) - cb.*sqrt(r.*p1.*(1-p1) + p0.*(1-p0)))./da;
n1 = s.*s./r;
n1 = n1.*(.5 + sqrt(.25 + (1+r)./(2*n1.*r.*da)))^2;
retp(n1);
endp;
