/* This proc performs a Wald test of each column of b[i,.] using
   the covariance matrix bcov[i,i].
   Inputs: b = hypothesis matrix, each column of which is b-hat minus
               a test value for b,
           bcov = covariance matrix for b[.,j] (assumed same for all j),
           names = names for each test (if scalar, same name will be used
                   for each test - set to zero for no printed output).
   Output: wstat = Wald statistic, one for each column of b,
           p = Wald P-values, one for each column of b.
*/
proc (2) = waldtest(b,bcov,names);
local df,wstat,p,j;
df = rows(b);
trap 1; wstat = diag(b'invpd(bcov)*b);
   if scalerr(wstat); df = rank(bcov);
      "SINGULARITY IN WALDTEST.G: df = "$+ftocv(df,1,0);
      wstat = diag(b'invpd(bcov)*b);
   endif; trap 0;
   if scalerr(wstat); retp(-1,-1); endif;
p = cdfchic(wstat,df);
if 0$+names[1]; j = 0;
  do until j ge cols(b); j = j + 1;
     names[j]$+": Wald ";; call chreport(wstat[j],df); endo;
endif;
retp(wstat,p);
endp;
