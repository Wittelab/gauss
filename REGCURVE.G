/* This proc computes and (if reqested) graphs the fitted values and their
   confidence limits for the given regressor values.
   Inputs: x = matrix of regressor values
               (fitted value will be computed for each row)
           abscis = column vector of values to serve graph abscissa
                  -- set to 0 for no graph -- requires name for graph
           const = 1 if constant to be added to x, 0 if not
           tabl = 1 to print table of fitted values, CL, & regressors, 0 not
           link = link function for model:
                  0 = identity link (linear model)
                  1 = log link (exponential model)
                  2 = logit link (logistic model)
           b = estimated coefficients
           bcov = estimated covariance matrix for b --
                note that if b contains intercepts and is from a logistic
                model for case-control data, the intercept variance estimates
                (on the diagonal of bcov) should be corrected as done in AFCC.G
                (see Prentice & Pyke, Biometrika 1979)
           xnames = names of columns of x (0 for no printed output)
           ynames = name of outcome (0 for no printed output)
           abname = name of abscissa (0 for no graph)
   Outputs: fit = matrix with fitted y, lower & upper limits in cols. 1-3
   Global (may be set by user from calling program; default is 0.95):
           _clevel = confidence level (expressed as a proportion)
*/
proc regcurve(x,abscis,const,tabl,link,b,bcov,xnames,yname,abname);
declare matrix _clevel = .95;
local z,eta,se,fit,i,isx;
if const; x = ones(rows(x),1)~x; endif;
/* standard normal percentile for CL: */ z = cinvnorm((1-_clevel)/2);
/* std err of fitted linear predictors: */ se = sqrt(sumr(x.*(x*bcov)));
/* fitted linear predictors & CL: */ eta = x*b + (0~(-z)~z).*se;
/* transform to fitted mean: */
if link eq 0; fit = eta;
   elseif link eq 1; fit = exp(eta);
   elseif link eq 2; fit = expit(eta);
   else; "ERROR IN REGCURVE.G: link must be 0, 1, or 2."; end;
endif;
if 0$+xnames[1]; /* tabulate fit and x: */
   if tabl;
      " Fitted "$+yname$+" and ";; format /rds 2,0; 100*_clevel;;
      "percent limits at values of regressors"; format 8,3; $xnames';;":";
      i = 0;
      do until i ge rows(x); i = i + 1;
         fit[i,.];; " for regressor values ";
         x[i,1+const:cols(x)];
      endo;
   endif;
   if 0$+abname; /* plot fit against designated regressor: */
      library pgraph; _pdate =""; fonts("microb");
      let _pstype = 8 10 13; let _psymsiz = 2.5 2 2;
      let _plctrl = 1 1 1;
      let _plwidth = 2 7 7; let _pltype = 6 2 2;
      xlabel(abname); ylabel(yname);
      /* sorted indices of x: */ isx = sortind(abscis);
      xy(abscis[isx],fit[isx,.]);
   endif;
endif;
retp(fit);
endp;
