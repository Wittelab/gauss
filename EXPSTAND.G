/* This proc computes standardized rate ratios and differences
   from exponential-regression output.
   Inputs: d = person-time at each covariate level of standard dist.
           b = coefficients
           bcov = covariance matrix of b
           x1 = index ("exposed") design matrix for numerator rate
           x0 = reference ("unexposed") design matrix for denominator
           const = 1 if b[1] is intercept and x contains no constant,
                   0 otherwise
           name = name of index level (set to zero for no printed output).
   Outputs: sr1,sr0 = index and reference standardized rates
            lsrr,srd = standardized ln(ratio) and difference estimates
            sr1se,sr0se = estimated standard errors for sr1 and sr0
            lsrrse,srdse = estimated standard errors for ln(rr) and rd
*/
proc (8) = expstand(d,b,bcov,x1,x0,const,name);
local sumd,a1,a0,sr1,sr0,dr1,dr0,srv1,srv0,
      scov,lsrr,srd,lsrd,lsrrse,lsrdse,srrse,srdse,lims;
if const; x1 = ones(rows(x1),1)~x1; x0 = ones(rows(x0),1)~x0; endif;
/* sum of standard: */ sumd = sumc(d);
/* index case expectations: */ a1 = exp(x1*b).*d;
/* index standardized rate: */ sr1 = sumc(a1)/sumd;
/* derivative of sr1 wrt b: */ dr1 = (x1'a1)/sumd;
/* sr1 variance: */ srv1 = dr1'bcov*dr1;
/* reference case expectations: */ a0 = exp(x0*b).*d;
/* reference standardized rate: */ sr0 = sumc(a0)/sumd;
/* derivative of sr0 wrt b: */ dr0 = (x0'a0)/sumd;
/* sr0 variance: */ srv0 = dr0'bcov*dr0;
/* cov(sr1,sr0): */ scov = dr1'bcov*dr0;
/* ln(ratio): */ lsrr = ln(sr1/sr0);
/* difference: */ srd = sr1 - sr0;
/* se(ln(srr)): */
   lsrrse = sqrt(srv1/(sr1^2) + srv0/(sr0^2) - 2*scov/(sr1*sr0));
/* se(srd): */ srdse = sqrt(srv1 + srv0 - 2*scov);
if 0$+name[1]; /* print results: */ format /rdn 12,6;
    "Standardized index and reference risks or rates:";; sr1;;sr0;
    "Standardized ratio and 95% limits from ln-transform:";
       format 12,4; exp(lsrr + lsrrse*(0~(-1.96)~1.96)); format 12,6;
    "Standardized difference and 95% limits, no transform:";
       format 12,6; srd+(0~(-1.96)~1.96)*srdse;
endif;
retp(sr1,sr0,sqrt(srv1),sqrt(srv0),lsrr,srd,lsrrse,srdse);
endp;
