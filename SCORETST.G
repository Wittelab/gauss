/* This proc computes the score test for columns of x with coefficients
   constrained to zero in a generalized linear model linear in the natural
   parameter, such as binomial-logistic and Poisson-exponential models.
   Inputs: x = full (augmented) design matrix INCLUDING CONSTANT IF PRESENT
           df = degrees of freedom for test
                (no. of cols. of x excluded from fit)
           e = residual vectors from constrained fit to means
               (e = obs. no. cases - no. expected under constrained model)
           w = weight matrix on linear scale from constrained fit
               (w = variance of obs. no. cases under constrained model)
           name = name for test (set to zero for no printed output).
   Outputs: sstat = score statistic
            p = p value from this statistic
*/
proc (2) = scoretst(x,df,e,w,name);
local s,sstat,p;
/* score vector: */ s = x'e;
/* score statistic: */ trap 1; sstat = s'invpd(moment(sqrt(w).*x,0))*s;
   if scalerr(sstat); df = df - cols(x) + rank(x);
      "SINGULARITY IN SCORETST.G: df = "$+ftocv(df,1,0);
      sstat = s'pinv(moment(sqrt(w).*x,0))*s;
   endif; trap 0;
   if scalerr(sstat); retp(-1,-1); endif;
if 0$+name; name$+": Score ";; p = chreport(sstat,df);
   else; p = cdfchic(sstat,df);
endif;
retp(sstat,p);
endp;
