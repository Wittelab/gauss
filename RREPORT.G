/* This proc generates report of results of RR regression.
   If rows(b) = rows(names)+1, program will assume b[1] is the intercept.
   NOTE: Returned P-values & confidence limits are from normal distribution
         unless df is negative, in which case a t distribution with -df
         degrees of freedom is used.
   Inputs: b = estimated coefficients
           se = estimated standard errors
           bnames = names of coefficents
           yname = name of dependent variable
           dev = residual deviance (negative if not given)
           rss = residual sum of squares (negative if not given)
           df = residual degrees of freedom (negative for t distribution)
           rr = 2 if relative-risk report and tests of fit desired,
                1 if only relative-risk report desired, 0 otherwise
   Globals (may be set by user from calling program):
           _clevel = confidence level expressed as a proportion
                     (default is .95)
           _dtime = do not print date if 0 (default is 1)
           _only = 0 if print full results for all coefficients in bnames, 
                   else vector of indices to print only RR, limits, P

*/
proc (0) = rreport(b,se,bnames,yname,dev,rss,df,rr);
declare matrix _clevel = .95;
declare matrix _dtime = 1;
declare matrix _only = 0; 
local c,p,np,i,z,prr; print;
if rows(bnames) eq rows(b)-1; /* add name for constant: */
   bnames = "constant"|bnames; c = 2;
   elseif rows(bnames) eq rows(b); c = 1;
   else; "ERROR: CANNOT MATCH b NAMES TO b IN RREPORT.G"; retp; endif;
if sumc(se .lt 0);
   "ERROR: SE CANNOT BE NEGATIVE -- SE vector: "; se'; retp; endif;
np = rows(b);
("REGRESSING "$+(yname')$+" ON "$+ftocv(np-c+1,1,0)$+" COVARIATE");;
if np-c; "S";; endif; if c ge 2; " PLUS CONSTANT";; endif;
if _dtime; " -- "$+datestr(0)$+" at "$+timestr(0);; endif; ".";
if df lt 0; p = (se .ne 0).*cdftt(b./(se + (se .eq 0)),-df);
   " Using a t distribution with ";; -df;; " degrees of freedom --";
    else; p = (se .ne 0).*cdfnt(b./(se + (se .eq 0)));
          " Using a normal distribution --";
endif;
i = 0;
if not _only[1] and rr ne 1;
   "Estimated betas, standard errors, & Wald P values:";
   do until i ge np; i = i + 1; 
      format /lds 8,0; $bnames[i];;
      format /rds 15,5; b[i]~se[i];; 
      format /ldn 9,6; "       ";;
      if p[i] ge .1; format 9,2; elseif p[i] ge .01; format 9,3;
         elseif p[i] ge .001; format 9,4; elseif p[i] ge .0001; format 9,5;
      endif; p[i];
   endo;
endif;
print;
if rr eq 1 or rr eq 2; /* print CLs: */
   if df lt 0; z = invt(1-(1-_clevel)/2,-df);
      else; z = invnorm(1-(1-_clevel)/2); endif;
      prr = exp(b + (0~-z~z).*se);
   "Estimated ratios per unit regressors,";;
   format 3,0; 100*_clevel;;"% Wald limits";;
   if _only[1]; ", and P-values";; endif; ":";
   /* Delete constant from RR output if c = 2: */ i = c - 1;
   do until i ge np; i = i + 1; 
      if _only[1]; if prodc(i-c+1 .ne _only); continue; endif; endif;
      format /lds 8,0; $bnames[i];;
      format /rds 15,5; prr[i,1 2];; ",";; 
      format 14,4; prr[i,3];;
      if _only[1]; format /ldn 9,6; "   ";;
         if p[i] ge .1; format 9,2; elseif p[i] ge .01; format 9,3;
            elseif p[i] ge .001; format 9,4; elseif p[i] ge .0001; format 9,5;
         endif; p[i];
         else; print;
      endif;
   endo;
endif; 

if (dev ge 0 or rss ge 0) and not _only[1]; 
   df = abs(df);
   if round(df) eq df and df ne 0;
         format /ldn 8,0; "Residual degrees of freedom : ";; df;
         elseif df ne 0;
         format /ldn 10,2; "Approximate residual degrees of freedom : ";; df;
    endif;
    if dev ge 0;
      format /ldn 10,2; "Final residual deviance: ";; dev;;
      if rr ge 2 and df gt 0;
         format /ldn 10,4; ", p = ";; cdfchic(dev,df);; 
      endif; print;
    endif;
    if rss ge 0;
      format /ldn 10,2; "Final weighted RSS     : ";; rss;;
      if rr ge 2 and df gt 0;
         format /ldn 10,4; ", p = ";; cdfchic(rss,df);; 
      endif; print;
    endif;
endif;
format /rds 9,3;
retp;
endp;
