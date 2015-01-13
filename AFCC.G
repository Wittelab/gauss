/* This procedure computes the attributable fraction and its variance from
   case-control data modeled by logistic regression, two ways:
   maximum likelihood (Greenland & Drescher, 1993); and Benichou & Gail, 1990.
   WARNING: ML METHOD ASSUMES THE FOLLOWING:
            1) FIRST ncons COEFFICIENTS ARE INTERCEPTS.
            2) INCOMING bcov MATRIX IS INVERSE-INFORMATION MATRIX,
               UNCORRECTED FOR CASE-CONTROL SAMPLING, as one would obtain
               from an ordinary logistic program (proc will correct intercept
               variances by subtracting 1/m1k + 1/m0k from them).
            3) ASSUMES UNMATCHED SAMPLING.
   Inputs: y = no. of cases at each exposure-covariate level
           n = no. of subjects at each exposure-covariate level
               (if scalar, all n will be set to this value)
           b = logistic coefficients
           bcov = covariance matrix of b
           x = observed design matrix for regression
           xr = reference design matrix
           ncons = number of constants in x (if cols(x) eq rows(b)-1,
                   proc will add constant to front of x; ncons should
                   then be 1)
           name = name of index level (set to zero for no printed output).
   Outputs: gaf,baf = ML & Bruzzi logistic attributable-fraction estimates
            gafse,bafse = estimated ML & B-G standard errors
*/
proc (4) = afcc(y,n,b,bcov,x,xr,ncons,name);
local c,k,m1,m0,d,r,s,p,gaf,baf,gafse,bafse,
      ps,psd,dn,dt,q,dp,dq,sdp,sdq,dc,dcx,dcxr;
n = n.*ones(rows(x),1);
/* Augment x with constant if x is 1 column less than rows(b): */
if cols(x) eq (rows(b)-1);
   x = ones(rows(x),1)~x; d = zeros(rows(x),1)~d;
   elseif cols(x) ne rows(b);
          "UNABLE TO MATCH COLS(x) TO ROWS(b) IN AFCC.G"; end;
endif;
/* save inverse-information: */ c = bcov;
/* correct intercept variances: */ k = 0;
   do until k ge ncons; k = k + 1;
      bcov[k,k] = bcov[k,k] -  1/sumc(y.*x[.,k]) - 1/sumc((n-y).*x[.,k]);
   endo;
/* count cases and noncases: */ m1 = sumc(y); m0 = sumc(n) - m1;
/* difference of reference and index levels: */ d = xr - x;
/* pseudo-risk estimates: */ r = expit(x*b);
/* inverse RR estimates: */ s = exp(d*b);
/* fitted distribution among cases: */ p = (n/m1).*r;
/* contributions to proportionate change: */ ps = p.*s;
/* ML estimate of AF: */ gaf = 1 - sumc(ps);
psd = d'ps;
dn = (r.*s)/m1;
dt = psd + x'((1-r).*ps);
q = (n/m0).*(1-r);
dp = dn.*p; dq = dn.*q; sdp = sumc(dp); sdq = sumc(dq);
dc = dt'c;
dcx = dc*x';
dcxr = dcx.*r';
/* Estimated variance of gaf: */
gafse = dt'bcov*dt
 + 2*m1*((dcx-dcxr)*dp - sumc((dcx-dcxr)'.*p)*sdp)
 - 2*m0*(dcxr*dq - sumc(dcxr'.*q)*sdq)
 + dn'(m1*dp + m0*dq) - m1*sdp*sdp - m0*sdq*sdq;
/* Estimated standard error of gaf: */ gafse = sqrt(gafse);
/* This code is clearer but runs much more slowly (~twice the time):
local va,vn,u;
va = m1*(diagm(p) - p*p');
vn = va + m0*(diagm(q) - q*q');
u = (va - vn.*r')*x*c;
gafse = sqrt(dt'bcov*dt + 2*(dn'u*dt) + dn'vn*dn);
*/
/* Bruzzi-Benichou-Gail estimates - note that the variances and covariances of
   intercept and confounder coefficients do not actually get used by the
   formula, because their columns of d are zero; however, their covariances
   with exposure coefficients do get used. */
local csp,psdc;
/* observed case distribution: */ p = y/m1;
/* estimated inverse of psmr: */ ps = p.*s;
/* Bruzzi estimate of AF: */ baf = 1 - sumc(ps);
psd = d'ps;
psdc = psd'c;
csp = (psdc)*(x'*((1-r).*(ps - p*sumc(ps))));
/* Estimated standard error of baf: */
   bafse = sqrt(psd'bcov*psd + 2*csp + (ps's - sumc(ps)^2)/m1);
/* Quadratic form of B&G variance:
   local v,cp; cp = diagm(p) - p*p';
   v = psd'bcov*psd + 2*m1*(psd'c*x'diagm(1-r)'cp*s) + m1*(s'cp*s);
   sqrt(v);
*/
format /rdn 10,4;
if 0$+name[1]; /* print results: */ local u,pc;
   "Estimated pop. attributable fraction AF & proportionate change PC = 1-AF";
   "  for "$+name$+",  based on logistic model from case-control data -";
   "By maximum-likelihood  : AF = ";; gaf;;" , PC = ";; pc = 1-gaf; pc;
   " 95% confidence limits for AF & proportionate change based on log(1-AF) -";
   u = exp(1.96*gafse/pc);
   1-(u~(1/u))*pc;;" ,";; ((1/u)~u)*pc;
   "By Bruzzi-Benichou-Gail: AF = ";; baf;;" , PC = ";; pc = 1-baf; pc;
   " 95% confidence limits for AF & proportionate change based on log(1-AF) -";
   u = exp(1.96*bafse/pc);
   1-(u~(1/u))*pc;;" ,";; ((1/u)~u)*pc;
endif;
retp(gaf,gafse,baf,bafse);
endp;
