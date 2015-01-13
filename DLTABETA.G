/* This proc does delta-beta analysis.
   Inputs: db = delta-betas, b = betas, se = standard errors
           x = design matrix,
           bnames = coefficient names (MUST be supplied).
           xp = covariates to plot delta-betas against;
                0 if no plot desired.
           xnames =  names of columns of xp; may be 0 if xp is 0.
*/
proc (0) = dltabeta(db,b,x,bnames,xp,xnames);
library pgraph;
external matrix _pnotify,_plctrl,_psymsiz,_pcrop;
/* put no curve through plots: */  _plctrl = -1;
/* use small plot symbols: */  _psymsiz = 2;
/* do not crop points outside axes: */ _pcrop = { 0 0 0 0 0 };
local se,i,m,j;
if rows(db) ne rows(b);
   "ERROR: CANNOT MATCH DELTABETAS TO b IN DLTABETA.G"; retp;
endif;
se = sqrt(diag(bcov));
/* 1. Remove first element of b if : */
if rows(bnames) eq rows(b)-1;
   b = b[2:rows(b)]; se = se[2:rows(se)]; db = db[.,2:cols(db)];
   elseif rows(bnames) ne rows(b);
          "ERROR: CANNOT MATCH b NAMES TO b ENTRIES IN DLTABETA.G"; retp;
endif;
/* 4. Summarize delta-betas: */
i = 0;
do until i ge cols(db); i = i + 1;
   format /rdn 10,4; "Summary of delta-betas for ";; $bnames[i];
   "   with beta and se of ";; b[i];; " and ";; se[i];;":";
   " Most positive delta-beta = ";; m = maxindc(db[.,i]); db[m,i];;
   " at covariate pattern";; format 6,0; m;; ":";
   format 10,4; $bnames'; x[m,.];
   " Most negative delta-beta = ";; m = minindc(db[.,i]); db[m,i];;
   " at covariate pattern";; format 6,0; m;; ":";
   format 10,4; $bnames'; x[m,.];
endo;
if rows(xp) gt 1; /* 3. Plot delta-betas: */
   j = 0;
   do until j ge cols(xp); j = j + 1;
      i = 0;
      do until i ge rows(b); i = i + 1;
         "Press enter for plot of delta-betas against ";;
         $xnames[j];;" for coefficient of ";; $bnames[i];
         wait; xy(xp[.,j],db[.,i]);
      endo;
  endo;
endif;
retp;
endp;
