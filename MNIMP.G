/* This proc imputes values for those missing in each column of a matrix xm
   by randomly sampling values from a fitted multivariate normal regression
   model for xm, using any complete variables xc as the regressors. To provide
   reasonably fast computing times, the missing-data means are fit by a
   noniterative least squares method rather than by maximum likelihood.
   Inputs: xc = variables with no missing values (enter 0 if there are none)
           WARNINGS: xc'xc must be invertible; do not include a constant in xc
           xm = variables with missing values
           mx = missing value indicator matrix for xm (mx has a 1 entered
                wherever xm has a missing value, and 0 entered elsewhere)
           nimp = number of imputations desired (ordinarily 2 or more)
                  -- if 0 is entered, the program will return 0 for ximp.
   Outputs: xm = xm with predicted values (conditional means) substituted for
                 missing values
            ximp = nimp stacked matrices, each identical to xm for
                   observed values and with random multivariate normal
                   imputations in place of missing values (0 if nimp = 0).
   NOTE: If ximp is to be used for further analyses, e.g. as input to the
         regression procedures LRIMP.G or CLRIMP.G, xc and xm should contain
         any variables potentially useful for imputing values of the variables
         of scientific interest, including regressands (variables that will
         later be regressed on columns of xm); for example, xc should include
         disease even if disease is to be regressed on  xm in later analyses
         (Little, JASA 1992, p. 1235). If however the returned xm is to be used
         for further analyses, regressands should be excluded from xc and xm.
*/
proc (2) = mnimp(xc,xm,mx,nimp);
local t0,nr,nmc,nm,c,m,b,xp,ruc,j,res,sd,rcov,ximp,s,ci,i,vu,vo,a,rcc;
t0 = date; "Proc mnimp.g: Multivariate normal multiple imputation.";
nr = rows(xm); nmc = cols(xm);
if rows(xc) eq nr; xc = ones(nr,1)~xc; else; xc = ones(nr,1); endif;
/* no. missing vars: */ nm = sumr(mx);
/* indices of complete subjects: */ c = selif(seqa(1,1,nr),nm .eq 0);
/* indices of incomplete subjects: */ m = selif(seqa(1,1,nr),nm .gt 0);
/* regress xm on xc using complete subjects: */
   /* compute coefficients for multivariate-outcome regression: */
      b = invpd(xc[c,.]'xc[c,.])*(xc[c,.]'xm[c,.]);
   /* unconditional predicted values: */ xp = xc*b;
/* estimate unconditional correlation matrix of residuals: */
   ruc = corrx(xm[c,.] - xp[c,.]);
/* estimate residual SDs from univariate-outcome regressions: */
   j = 0; sd = zeros(nmc,1);
   do until j ge nmc; j = j + 1;
      c = selif(seqa(1,1,nr),1-mx[.,j]);
      res = xm[c,j] - xc[c,.]*invpd(xc[c,.]'xc[c,.])*(xc[c,.]'xm[c,j]);
      sd[j]  = sqrt((res'res)/(rows(c)-rows(b)));
   endo;
   /* unconditional covariance matrix: */ rcov = sd.*ruc.*sd';
/* do imputations based on the fitted multinormal model: */
   if nimp;  /* initialize ximp to nimp copies of xm: */
      ximp = xm[ones(nimp,1).*.seqa(1,1,nr),.];
      /* root unconditional covariance: */ ruc = chol(rcov);
      /* start points for indices in ximp: */ s = seqa(0,nr,nimp);
      else; ximp = 0; ruc = 0;
   endif;
   j = 0; ci = seqa(1,1,nmc);
   do until j ge rows(m); j = j + 1; i = m[j];
      /* impute missing values for subject i by multinormal regression: */
      /* indices of unobserved vars: */ vu = selif(ci,mx[i,.]');
      if nm[i] lt nmc; /* regress on observed values in xm[i,.]: */
         /* indices of observed vars: */ vo = selif(ci,1-mx[i,.]');
         a = invpd(rcov[vo,vo])*rcov[vo,vu];
         /* substitute conditional means for missing values: */
            xm[i,vu] = xp[i,vu] + (xm[i,vo] - xp[i,vo])*a;
         if nimp; /* compute root conditional covariance: */
            rcc = chol(rcov[vu,vu] - rcov[vu,vo]*a);
         endif;
         else; /* no obs values in xm[i,.], so use unconditional dist; */
               xm[i,.] = xp[i,.]; rcc = ruc;
      endif;
      if nimp; /* generate random values: */
         ximp[s+i,vu] = xm[i,vu] + rndn(nimp,rows(vu))*rcc; endif;
   endo;
"Total imputation time: ";; etstr(ethsec(t0,date));
retp(xm,ximp);
endp;
