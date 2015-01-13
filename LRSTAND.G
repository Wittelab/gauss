/* This proc computes standardized proportion ratios and differences
   from logistic-regression output.
   Inputs: d = no. or proportion at each covariate level of standard dist.
           b = logistic coefficients
           bcov = covariance matrix of b
           x1 = index ("exposed") design matrix for numerator proportion
           x0 = reference ("unexposed") design matrix for denominator
           const = 1 if b[1] is intercept and x contains no constant,
                   0 otherwise
           name = name of index level (set to zero for no printed output).
   Outputs: sp1,sp0 = index and reference standardized proportions
            lspr,spd = standardized ln(ratio) and difference estimates
            sp1se,sp0se = estimated standard errors for sp1 and sp0
            lsprse,spdse = estimated standard errors for ln(rr) and rd
*/
proc (8) = lrstand(d,b,bcov,x1,x0,const,name);
local sumd,p1,p0,a1,a0,sp1,sp0,dp1,dp0,spv1,spv0,
      scov,lspr,spd,lspd,lsprse,lspdse,sprse,spdse,lims;
if const; x1 = ones(rows(x1),1)~x1; x0 = ones(rows(x0),1)~x0; endif;
/* sum of standard: */ sumd = sumc(d);
/* index p: */ p1 = expit(x1*b);
/* index case expectations: */ a1 = p1.*d;
/* index standardized proportion: */ sp1 = sumc(a1)/sumd;
/* derivative of sp1 wrt b: */ dp1 = (x1'(a1.*(1-p1)))/sumd;
/* sp1 variance: */ spv1 = dp1'bcov*dp1;
/* reference p: */ p0 = expit(x0*b);
/* reference case expectations: */ a0 = p0.*d;
/* reference standardized proportion: */ sp0 = sumc(a0)/sumd;
/* derivative of sp0 wrt b: */ dp0 = (x0'(a0.*(1-p0)))/sumd;
/* sp0 variance: */ spv0 = dp0'bcov*dp0;
/* cov(sp1,sp0): */ scov = dp1'bcov*dp0;
/* ln(ratio): */ lspr = ln(sp1/sp0);
/* difference: */ spd = sp1 - sp0;
/* se(ln(spr)): */
   lsprse = sqrt(spv1/(sp1^2) + spv0/(sp0^2) - 2*scov/(sp1*sp0));
/* se(spd): */ spdse = sqrt(spv1 + spv0 - 2*scov);
format /rdn 12,6;
if 0$+name[1]; /* print results: */
    "Standardized index and reference proportions:";; sp1;;sp0;
    "Standardized ratio and 95% limits from ln-transform:";
       format 12,4; exp(lspr);; exp(lspr + 1.96*lsprse*(-1~1)); format 12,6;
    "Standardized difference and 95% limits from logit transform:";
       /* logit transform for spd: */ lspd = ln((1 + spd)/(1 - spd));
       /* se(lspd): */ lspdse = 2*spdse/((1 + spd)*(1 - spd));
       lims = lspd + 1.96*lspdse*(-1~1);
       /* transform limits back to spd scale: */
          spd;; (exp(lims) - 1)./(exp(lims) + 1);
endif;
retp(sp1,sp0,sqrt(spv1),sqrt(spv0),lspr,spd,lsprse,spdse);
endp;
