/* This proc computes covariate-specific risk or rate ratios and differences
   from exponential (log-linear) regression output.
   Inputs: b = log-linear coefficients
           bcov = covariance matrix of b
           x1 = index ("exposed") design matrix for numerator risk or rate
           x0 = reference ("unexposed") design matrix for denominator
           const = 1 if b[1] is intercept and x contains no constant,
                   0 otherwise
           name = name of index level (set to zero for no printed output).
   Outputs: r1,r0 = index and reference specific risks or rates
            lrr,rd = specific ln(ratio) and difference estimates
            vr1,vr0 = estimated covariance matrices for p1 and p0
            vlrr,vrd = estimated covariance matrices for ln(rr) and rd
*/
proc (8) = expspec(b,bcov,x1,x0,const,name);
local r1,r0,dr1,dr0,vr1,vr0,cov,lrr,rd,vlrr,vrd,rlims,dlims;
if const; x1 = ones(rows(x1),1)~x1; x0 = ones(rows(x0),1)~x0; endif;
/* index risks or rates: */ r1 = exp(x1*b);
/* derivative of r1 wrt b: */ dr1 = x1.*r1;
/* r1 covariance matrix: */ vr1 = dr1*bcov*dr1';
/* reference risks or rates: */ r0 = exp(x0*b);
/* derivative of r0 wrt b: */ dr0 = x0.*r0;
/* r0 covariance matrix: */ vr0 = dr0*bcov*dr0';
/* cov(r1,r0): */ cov = dr1*bcov*dr0';
/* ln(ratio): */ lrr = ln(r1./r0);
/* difference: */ rd = r1 - r0;
/* cov(ln(rr)): */
   vlrr = vr1./(r1^2) + vr0./(r0^2) - 2*cov./(r1.*r0);
/* cov(vrd): */ vrd = vr1 + vr0 - 2*cov;
if 0$+name[1]; /* print results: */ "PROC EXPSPEC.G:";;
   " Exposure-specific risks or rates from exponential regression.";
   "Fitted index and reference risks or rates,";
   "ratios and 95% limits from ln-transform,";;
   "and differences and 95% limits:"; format 9,4;
    rlims = exp(lrr + sqrt(diag(vlrr))*(0~(-1.96)~1.96));
    dlims = rd + sqrt(diag(vrd))*(0~(-1.96)~1.96); r1~r0~rlims~dlims;
endif;
retp(r1,r0,vr1,vr0,lrr,rd,vlrr,vrd);
endp;
