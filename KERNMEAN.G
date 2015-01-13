/* This proc computes kernel (running weighted) means using the Epanechnikov
   kernel function.
Inputs: xf = matrix of points at which means will be computed
             -- set to 0 if xf = x
        h = kernel radii (scalar, or matrix or vector conformable with xf)
        x = design matrix (regressor matrix)
            -- it is assumed that data have been grouped so that
               the rows of x are UNIQUE, and that the distances between
               rows of x are MEANINGFUL (so that the columns of x are
               quantitative, not category indicators)
        y = outcome vector (regressand)
        n = denominators for Poisson or binomial y, repetition count otherwise
            (may be a scalar if all counts are equal)
        v = vector of variances for y,
            OR code for variance function for y:
            v = 1 if y is Poisson, in which case v = yhat will be used
                2 if y is binomial, in which case v = yhat.*(1 - yhat./n)
        bnames = names for columns of x (set to zero for no printed output)
        yname = name for y
Output: yf = fitted values for y./n (if v = 1 or 2) or for y
        yfcov = estimated covariance matrix for yf (large-cell unless xf = x)
        dev = residual deviance if y is Poisson or binomial
        rss = residual sum of squares, weighted by 1/v
*/
proc (4) = kernmean(xf,h,x,y,n,v,bnames,yname);
local gfit,nf,k,kn0,yf,veq0,i,r,e,dev,rss,nn0,nz;
gfit = rows(xf) eq 1 and cols(xf) eq 1 and xf[1,1] eq 0;
if gfit; xf = x; endif;
nf = rows(xf);
if rows(n) eq 1; n = n*ones(rows(x),1); endif;
if rows(h) eq 1; h = h[ones(nf,1),.]; endif;
/* construct weight matrix: */
k = zeros(nf,rows(x)); i = 0;
do until i ge nf; i = i+1; /* squared distances from point i: */
   k[i,.] = sumr(((x-xf[i,.])./h[i,.])^2)';
endo;
k = max(cols(x) - k,0); kn0 = k*n .le 0; k = (k.*n')./(k*n + kn0);
dev = 0; rss = 0; nn0 = n + (n .eq 0);
if rows(v) eq 1; r = k*(y./nn0); yf = r;
   nz = sumc(yf .eq 0 .or yf .eq 1 .or kn0);
   if nz; "WARNING:";; format 4,0; nz;;" fitted values equal 0 or 1 --";;
          " increase radii h."; format 9,3;
   endif;
   else; yf = k*y;
endif;
if gfit;
   if rows(v) eq 1; e = n.*r;
      if v eq 1; v = e; veq0 = v .eq 0;
         dev = -2*sumc(-e + y.*ln(e+veq0) - lnfact(y+veq0));
         elseif v eq 2; v = e.*(1-r); veq0 = v .eq 0;
                dev = -2*sumc(y.*ln(r+veq0) + (n-y).*ln(1-r+veq0)
                      + lnfact(n) - lnfact(y) - lnfact(n-y));
      endif;
      rss = ((1-veq0).*(y-e))'((y-e)./(v+veq0)); v = v./nn0^2;
      else; e = yf; veq0 = v .eq 0;
            rss = ((1-veq0).*(y-e))'((y-e)./(v+(v.eq 0)));
   endif;
   if 0$+bnames[1]; "Deviance and RSS: ";; format 10,3; dev;; rss; endif;
   elseif rows(v) eq 1;
      if v eq 1; v = y./nn0^2;
         elseif v eq 2; v = y.*(n-y)./nn0^3; endif;
endif;
retp(yf,k*(v.*k'),dev,rss);
endp;
