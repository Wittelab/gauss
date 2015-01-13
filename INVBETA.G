/* This proc computes quantiles of the beta(a,b) distribution
   using Newton's method on the logit scale.
Inputs: q = desired quantiles (percentiles/100)
        a,b = beta parameters
Return: p = the points such that cdfbeta(p,a,b) = q
   Globals: _maxit = max. no. of iterations (default: 20)
            _criter = convergence criterion (default: 0.0025)
*/
proc invbeta(q,a,b);
declare matrix _maxit = 20; declare matrix _criter = .0025;
local x,xold,i,cnv;
if sumall(q .lt 0 .or q .gt 1);
   "INVBETA INPUT ERROR: q < 0 OR q > 1"; end; endif;
/* intialize by logistic-normal approx. : */
   x = ln(a./b) + invnorm(q).*sqrt(1/a + 1/b);
i = 0; cnv = 0;
do until i ge _maxit or cnv; i = i + 1; xold = x;
   x = x - (cdfbeta(expit(x),a,b) - q)./lgstbeta(x,a,b);
/* q;; i;; xold;; x;;
cdfbeta(expit(x),a,b);; lgstbeta(x,a,b); */
   cnv = prodall(abs(exp(xold-x)-1) .lt _criter);
endo;
if not cnv; "WARNING: NO CONVERGENCE IN INVBETA.G AFTER";;
            _maxit;; " ITERATIONS."; endif;
retp(expit(x));
endp;
