/* This proc computes covariate-specific proportion ratios and differences
   from logistic-regression output.
   Inputs: b = logistic coefficients
           bcov = covariance matrix of b
           x1 = index ("exposed") design matrix for numerator proportion
           x0 = reference ("unexposed") design matrix for denominator
           const = 1 if b[1] is intercept and x contains no constant,
                   0 otherwise
           name = name of index level (set to zero for no printed output).
   Outputs: p1,p0 = index and reference specific proportions
            lpr,pd = specific ln(ratio) and difference estimates
            vp1,vp0 = estimated covariance matrices for p1 and p0
            vlpr,vpd = estimated covariance matrices for ln(rr) and rd
*/
proc (8) = lrspec(b,bcov,x1,x0,const,name);
local p1,p0,dp1,dp0,vp1,vp0,cov,lpr,pd,vlpr,vpd;
if const; x1 = ones(rows(x1),1)~x1; x0 = ones(rows(x0),1)~x0; endif;
/* index p: */ p1 = expit(x1*b);
/* derivative of p1 wrt b: */ dp1 = x1.*(p1.*(1-p1));
/* p1 covariance matrix: */ vp1 = dp1*bcov*dp1';
/* reference p: */ p0 = expit(x0*b);
/* derivative of p0 wrt b: */ dp0 = x0.*(p0.*(1-p0));
/* p0 covariance matrix: */ vp0 = dp0*bcov*dp0';
/* cov(p1,p0): */ cov = dp1*bcov*dp0';
/* ln(ratio): */ lpr = ln(p1./p0);
/* difference: */ pd = p1 - p0;
/* cov(ln(pr)): */
   vlpr = vp1./(p1^2) + vp0./(p0^2) - 2*cov./(p1.*p0);
/* cov(vpd): */ vpd = vp1 + vp0 - 2*cov;
if 0$+name[1]; /* print results: */ local rlims,lpdse,lpd,dlims;
    "PROC LRSPEC.G: Exposure-specific proportions from logistic regression.";
    "Fitted index and reference proportions,";
    "ratios and 95% limits from ln-transform,";;
    "and differences and 95% limits from logit transform:"; format 9,4;
       rlims = exp(lpr + sqrt(diag(vlpr))*(0~(-1.96)~1.96));
       /* logit transform for difference: */ lpd = ln((1 + pd)./(1 - pd));
       /* se(lpd): */ lpdse = 2*sqrt(diag(vpd))/((1 + pd).*(1 - pd));
       dlims = lpd + lpdse*(0~(-1.96)~1.96);
       /* transform limits back to difference scale: */
          dlims = (exp(dlims) - 1)./(exp(dlims) + 1);
       p1~p0~rlims~dlims;
endif;
retp(p1,p0,vp1,vp0,lpr,pd,vlpr,vpd);
endp;
