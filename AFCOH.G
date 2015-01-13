/* This procedure computes the attributable fraction and its variance from
   cohort data modeled by logistic regression, in two ways:
   iafse is based on score-beta covariance formula from Cox & Hinkley,
         cov(mx,b) = (1,0,...,0) if intercept is first term in b,
         where mx is observed total no. cases.
   See Greenland & Drescher, Biometrics, 1983.
   WARNING: iafse assumes first column in x is constant!
   gafse does not assume this, and does not use score-beta covariance formula.
   Inputs: y = no. of cases at each exposure-covariate level
           n = no. of subjects at each exposure-covariate level
           b = logistic coefficients
           bcov = covariance matrix of b
           x = observed design matrix for regression
           xr = reference design matrix
           const = 1 if b[1] = constant and x contains no constant,
                   0 if b[1] eq intercept and first column of x is constant
           name = name of index level (set to zero for no printed output).
   Outputs: af = logistic attributable-fraction estimate
            iafse,gafse = estimates of standard error
*/
proc (3) = afcoh(y,n,b,bcov,x,xr,const,name);
local rz,az,mx,mz,wz,af,pc,rx,wx,uz,iafse,ux,gafse,u;
if const ne 0; x = ones(rows(x),1)~x; xr = ones(rows(xr),1)~xr; endif;
/* fitted exposure-covariate distribution among cases: */
/* total no. cases: */ mx = sumc(y);
/* reference risks: */ rz = expit(xr*b);
/* reference cases expectations: */ az = rz.*n;
/* reference total expected: */ mz = sumc(az);
/* estimated proportionate change and AF: */ pc = mz/mx; af = 1 - pc;
/* reference case variances: */ wz = az.*(1 - rz);
/* index risks: */ rx = expit(x*b);
/* index case variances: */ wx = n.*rx.*(1 - rx);
/* Two types of standard error: */
uz = xr'wz;
iafse = sqrt(uz'bcov*uz + mz*(-2*sumc(wz) + mz*sumc(wx)/mx)/mx)/mx;
ux = (xr|x)'(wz|(-wx*mz/mx));
gafse = sqrt(ux'bcov*ux)/mx;
format /rdn 10,4;
if 0$+name[1]; /* print results: */
   "Estimated cohort attributable fraction AF & proportionate change PC = 1-AF";
   "  for "$+name$+",  based on logistic model from cohort data:";
   "AF = ";; af;;" , PC = ";; pc;
   "95% confidence limits for AF & proportionate change based on log(1-AF) -";
   "    from general formula        :";;
   u = exp(1.96*gafse/pc);
   1-(u~(1/u))*pc;;" ,";; ((1/u)~u)*pc;
   "    from intercept-first formula:";;
   u = exp(1.96*iafse/pc);
   1-(u~(1/u))*pc;;" ,";; ((1/u)~u)*pc;
endif;
retp(af,gafse,iafse);
endp;
