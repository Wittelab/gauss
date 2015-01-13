/* This proc does residual and deltabeta plots.
   Inputs: eta = fitted predictor vector.
           res = residual vector.
           sres = studentized residual vector.
           xp = matrix of covariates to plot residuals & deltabetas against;
                0 if no covariate-residual plots desired.
           db = deltabetas; 0 if no deltabeta plots desired.
           xnames = names of columns of xp; may be 0 if xp is 0.
           bnames = names of columns of db; may be 0 if db is zero.
*/
proc (0) = diagplot(eta,res,sres,xp,db,xnames,bnames);
library pgraph; graphset;
external matrix _pnotify,_plctrl,_psymsiz._pcrop;
/* put no curve through plots: */  _plctrl = -1;
/* use small plot symbols: */  _psymsiz = 2;
/* do not crop points outside axes: */ _pcrop = { 0 0 0 0 0 };
local i,j;
/* 1. Do residual plots: */
    "Press enter for histogram of studentized residuals.";
    "Ordinarily, this should have standard normal shape and spread.";
    xlabel("Studentized residuals"); wait; call hist(sres,10);
    "Press enter for plot of studentized residuals against linear predictor.";
    "Wrong choice of variance function may produce a trend in spread.";
    "Wrong choices of link or covariate forms may produce patterns";
    xlabel("Linear predictor"); ylabel("Studentized residuals");
    wait; xy(eta,sres);
    "Press enter for plot of linear-scale residuals against linear predictor.";
    "Wrong choices of link or covariate forms may produce patterns";
    xlabel("Linear predictor"); ylabel("Residuals"); wait; xy(eta,res);
if rows(xp) gt 1; /* 2. Do residual plots against covariates in xp: */
   if cols(xp) ne rows(xnames);
      "ERROR: CANNOT MATCH xp COLUMNS TO xnames IN DIAGPLOT.G"; retp; endif;
   if cols(db) ne rows(bnames);
      "ERROR: CANNOT MATCH db COLUMNS TO bnames IN DIAGPLOT.G"; db = 0; endif;
   j = 0;
   do until j ge cols(xp); j = j + 1;
      "Press enter for plot of residuals against ";; $xnames[j];
      "Wrong choice of covariate form will produce deviations from horizontal.";
      xlabel(xnames[j]); ylabel("Residuals");
      wait; xy(xp[.,j],res); ylabel("Delta-betas");
      if rows(db) gt 1; /* 3. Do deltabeta plots against covariates in xp: */
         i = 0;
         do until i ge cols(db); i = i + 1;
            "Press enter for plot of delta-betas against ";;
            $xnames[j];;" for coefficient of ";; $bnames[i];
            wait; xy(xp[.,j],db[.,i]);
         endo;
      endif;
   endo;
endif;
retp;
endp;

