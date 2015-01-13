/* This proc computes the conditional logistic deviance, score,
   and information at b for matched-set case-control sum differences xs.
   Inputs: xs = matched-set case-control covariate sum differences
           ofs = matched-set case-control offset sum differences
                 (0 if no offsets)
           NOTE: set-members should be contiguous in xs and ofs
           b = coefficents to be tested
           m = no. of differences for each matched set
               NOTE: sumc(m) should equal rows(xs)
           reps = repetion count for sets (set to vector of ones if none)
           ids = vector of matched-set indices in xs
   Note: cp = fitted probablity that the observed cases are the cases of the
              matched set, given the covariate and case totals.
   Outputs: dev = conditional deviance for b
            score = score vector at b
            inf = information matrix at b
*/
proc (3) = mscore(xs,ofs,b,m,reps,ids);
local r,s,cp,sm,dev,score,inf;
/* products of relative risks for each combination in each matched set: */
   r = exp(xs*b + ofs);
/* conditional probability of observed case pattern: */ 
   cp = 1/(1 + mpsumc(r,m));
/* scores for each combination: */ s = -cp[ids].*(xs.*r);
/* scores for each matched set: */ sm = mpsumc(s,m);
dev = -2*(reps'ln(cp));
score = sumc(reps.*sm);
inf = -(xs'(reps[ids].*s))-(reps.*sm)'sm;
retp(dev,score,inf);
endp;
