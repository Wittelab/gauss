/* This procedure does weighted-least squares regression of y on x, with an
   information-theory covariance estimate (w assumed to be inverse covariance
   matrix of y columns, or vector of inverse variances of y components).
   NOTE: OLS estimates can be obtained using w = 1, in which case each bcov
         will be premultiplied by its rss/df, the residual mean square, and
         a t distribution will be used for p-values & confidence limits.
   Inputs: x = design matrix,
           y = outcome vector or matrix (if matrix, b will be matrix),
           w = weight matrix or column vector (use scalar 1 for OLS),
           rep = repetition count (0 if none),
           const = 1 if constant term to be added to x, 0 otherwise
                   (if const = 0, a no-intercept model will be fitted
                    unless x already contains a constant column),
           bnames = names of regressors (set to 0 for no printed output),
           ynames = names of regressands (set to 0 for no printed output).
   Outputs: b = coefficient estimates,
            bcov = estimated covariance matrix for b
                   (NOT premultiplied by rss/rdf unless w = 1 or _res = 1,
                    in which case, if matrix b, covariance matrices of the
                    columns are vertically stacked
            rss = weighted residual sum of squares,
            r2 = weighted R-squared for regression, defined as one minus the
                  ratio of weighted RSS for the fitted model over the weighted
                  RSS for the null (no-covariate) model.
            rdf = residual degrees of freedom.
   Globals: _rdf = user-specified residual df; default n-cols(x).
           _res = 1 for bcov = (rss/df)*(x'w*x), 0 for bcov = x'w*x;
                  default 1 if w = 1, 0 otherwise.
           _rr = 1 for ratio confidence limits, 0 otherwise; default 0.
*/
proc (5) = wls(x,y,w,rep,const,bnames,ynames);
declare matrix _rdf, _res,_rr;
local ol,n,rdf,xw,b,bcov,res,yc,rss,r2,i,se;
/* indicate ols: */ ol = rows(w) eq 1;
if 0$+bnames[1]; print;
  if ol; "Proc WLS.G: Ordinary least squares regression.";
     else; "Proc WLS.G: Inverse-variance WLS regression";;
     if _res; " with residual variance."; else; "."; endif;
  endif;
  if rows(bnames) ne cols(x) or rows(ynames) ne cols(y);
     "INPUT ERROR: No. names supplied not equal to no. of variables."; end;
  endif;
endif;
if const; /* add constant: */ x = ones(rows(x),1)~x;
   /* center y for R-squared: */ yc = y - sumc(w'y)'/sumall(w); else; yc = y;
endif;
if rank(x'x) lt cols(x); "DESIGN MATRIX RANK DEFICIENT IN WLS.G";
   retp(0,0,0,0,0); endif;
if rows(rep) eq rows(x); /* use repetition counts: */
             w = rep.*w; n = sumc(rep); else; n = rows(x); endif;
if cols(w) eq 1; xw = (w.*x)'; else; xw = (w*x)'; endif;
/* inverse information: */ trap 1; bcov = invpd(xw*x); trap 0;
   if scalerr bcov; "SINGULARITY IN WLS.G"; retp(0,0,0,0,0); endif;
/* LS estimates: */ b = bcov*xw*y;
/* residual: */ res = y - x*b;
/* residual degrees of freedom: */  rdf = n - cols(x); if _rdf; rdf = _rdf; endif;
/* RSS and R-squared: */
   if cols(w) eq 1; /* weight vector: */
      rss = w'(res^2); r2 = 1 - rss./(w'(yc^2));
      else; /* weight matrix: */
            rss = diag(res'(w*res)); r2 = 1 - rss./diag(yc'(w*yc));
   endif;
if ol; /* premultiply each bcov by RMS: */ bcov = (rss/rdf).*.bcov; endif;
if 0$+bnames[1];
  if ol or _res;
     /* use negative df to indicate t distribution: */ rdf = -rdf; endif;
  format /rds 10,0; "No. observations:";; n;;
  i = 0;
  do until i ge cols(y); i = i + 1;
     se = sqrt(diag(bcov[(i-1)*rows(b)+seqa(1,1,rows(b)),.]));
     rreport(b[.,i],se,bnames,ynames[i],-1,rss[i],rdf,_rr*2); format 5,3;
     if not ol; "Weighted ";; endif; "R-squared: ";;  r2[i];;
     if not const; " (warning: no intercept added by program)."; else; print; endif;
  endo; format /rds 9,3;
endif;
retp(b,bcov,rss,r2,abs(rdf));
endp;
