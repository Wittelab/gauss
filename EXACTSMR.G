/* Test program: format 10,6; call exactsmr(10,1,"SMR"); end; */
/* This proc computes Fisher-exact confidence limits for an SMR or rate
   using the exact relation of the Poisson and chi-squared distributions.
   Inputs: a = no. observed cases
           e = no. cases expected OR person-time
           name = name for a/e (set to zero for no output)
   Outputs: ll,ul = lower and upper Fisher limits
   Globals: _clevel = confidence level (default is .95)
            _criter = convergence criterion (default is .0001)
            _maxit = maximum number of iterations (default is 30)
*/
proc (2) = exactsmr(a,e,name);
declare matrix _clevel = .95;
declare matrix _criter = .0001; declare matrix _maxit = 30;
local ha,z,ll,ul,lold,pold,cnv,iter,p;
ha = _clevel + (1-_clevel)/2; z = invnorm(ha);
/* initialize: */ ll = a*(1 - 1/(9*a) - z/(3*sqrt(a)))^3;
lold = (sqrt(a) - z/2)^2; pold = 1-cdfchic(2*lold,2*a);
cnv = 0; iter = 0;
do until cnv or iter ge _maxit; iter = iter+1;
   p = 1-cdfchic(2*ll,2*a);
   if abs(p - pold) lt _criter; cnv = 1;
      else; ll = ll - (ll-lold)*(p-1+ha)/(p-pold);
            lold = ll; pold = p;
   endif;
endo;
if iter gt _maxit;
   "WARNING: NO CONVERGENCE AFTER "$+ftocv(iter,2,0)$+"ITERATIONS."; endif;
/* initialize: */ ul = (a+1)*(1 - 1/(9*a+9) + z/(3*sqrt(a+1)))^3;
lold = (sqrt(a+1) + z/2)^2; pold = cdfchic(2*lold,2*a+2);
cnv = 0; iter = 0;
do until cnv or iter ge _maxit; iter = iter+1;
   p = cdfchic(2*ul,2*a+2);
   if abs(p - pold) lt _criter; cnv = 1;
      else; ul = ul - (ul-lold)*(p-1+ha)/(p-pold);
            lold = ul; pold = p;
   endif;
endo;
if 0$+name; format /rds 9,4;
   ("PROC EXACTSMR.G: Fisher-exact "$+ftocv(100*_clevel,2,0)$+
    " percent confidence limits");
    "                 for the "$+name$+":";; ll/e;;ul/e;
endif;
retp(ll/e,ul/e);
endp;
