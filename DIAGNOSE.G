/* This proc does studentized residual analysis for generalized-linear
   regessions with canonical links.
   Inputs: x = regressor (design) matrix, including constant if present
           y = observed outcome
           mu = mean function (must  be preceded by ampersand)
           wz = weight (variance of y) function (must be preceded by ampersand)
           b = coefficients
           bcov = estimated covariance matrix for b
           offset = vector or matrix of offsets (0 if no offsets)
           bnames = coefficient names; if 0, no output will be printed.
   Outputs: res = residuals on linear (z) scale
            sres = studentized residuals
            db = approximate (one-step) delta-betas (excluding intercept)
                 (db[i,j] is change in b[j]  from dropping row i of x)
 */
proc (4) = diagnose(x,y,&mu,&wz,b,bcov,offset,bnames);
local mu:fn,wz:fn,c,eta,muf,wf,e,w,
      np,cbsdx,sdy,ch,res,sres,db,z,m,i,j,se;
/* 1. If the no. of coefficients is one more than the number of columns of x,
   proc assumes first element of b is the intercept and so adds a column of
   ones to the front of x: */
   if cols(x) eq rows(b); c = 1;
      elseif (cols(x) + 1) eq rows(b); c = 2; x = ones(rows(x),1)~x;
      else; "ERROR: NO. COEFFICIENTS DOESN'T EQUAL COLS(x) OR COLS(x)+1"; end;
   endif;
np = cols(x);
/* 2. Compute fitted values and residuals: */
    /* linear predictor: */ eta = sumr(offset) + x*b;
    /* fitted mean: */ muf = mu(eta);
    /* fitted variance of y: */ wf = wz(eta);
    /* residual: */ e = y - muf;
    /* Standard deviation of y: */ sdy = sqrt(wf);
    /* h = sqrt(wf).*x*bcov*(x.*sqrt(wf))' may be too big to fit in memory, so
       use the following trick to avoid computing h. h = cbsdx'cbsdx but we
       don't need all of h, just its diagonal, which is sumc(cbsdx.*cbsdx): */
       cbsdx = chol(bcov)*(sdy.*x)';
    /* (CHOL takes the Cholesky square-root of a positive-definite matrix v,
       which is the upper-triangular matrix t such that t't = v). */
    /* variance adjustment factor for estimation of p: */
        ch = 1-sumc(cbsdx.*cbsdx);
    /* pseudo-residual in linear (z) scale for canonical links: */
        res = e./(wf + (wf .eq 0));
    /* studentized residual: */ sres = e./((sdy + (sdy .eq 0)).*sqrt(ch));
    /* approximate delta-betas: */ db = -(e./ch).*(x*bcov);
    /* exclude deltabeta for intercept: */ db = db[.,c:np];
if bnames[1];
   if rows(bnames) ne rows(b[c:np]);
      "ERROR: CANNOT MATCH b NAMES TO b ENTRIES IN DIAGNOSE.G";
   else; /* 3. Print residual summaries: */
      "Summary of studentized residuals:";
      " Most positive studentized residual =";; m = maxindc(sres);
      format /rdn 6,2; z = sres[m]; z; format 5,0;
      " at covariate pattern ";; m;; ":";
      format 10,4; $bnames'; x[m,c:np];
      format 10,7; "  Normal probability: ";; 1-cdfn(z);
      " Most negative studentized residual =";; m = minindc(sres);
      format 6,2; z = sres[m]; z; format 5,0;
      " at covariate pattern";; m;; ":";
      format 10,4; $bnames'; x[m,c:np];
      format 10,7; "  Normal probability: ";; cdfn(z);
   /* 4. Summarize delta-betas: */
      "Correlation matrix of linear-scale residuals and delta-betas -";
       format 10,3; "  residual";; if c eq 2; "  constant";; endif;
       $bnames'; corrx(res~db);
       se = sqrt(diag(bcov));
       i = 0;
       do until i ge cols(db); i = i + 1; j = i + c - 1;
          format /rdn 10,4; "Summary of delta-betas for ";; $bnames[i];
          "   with beta and se of ";; b[j];; " and ";; se[j];;":";
          " Most positive delta-beta = ";; m = maxindc(db[.,i]); db[m,i];;
          " at covariate pattern";; format 6,0; m;; ":";
          format 10,4; $bnames'; x[m,c:np];
          " Most negative delta-beta = ";; m = minindc(db[.,i]); db[m,i];;
          " at covariate pattern";; format 6,0; m;; ":";
          format 10,4; $bnames'; x[m,c:np];
       endo;
  endif;
endif;
/* 5. Return linear predictors & residuals, studentized residuals,
   and approximate delta-betas: */
retp(eta,res,sres,db);
endp;
