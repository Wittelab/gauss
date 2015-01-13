/* This procedure does conditional logistic regression on multiple-imputed
   1:m matched case-control data, using Rubin's method.
   A stack of multiple-imputed covariate matrices suitable for input to this
   procedure may be generated from MNIMP.G; however --
         MAKE SURE THAT THE ORIGINAL DATA ARE SORTED ON id,
         THE MATCHED-SET ID NUMBERS, BEFORE GENERATING THE IMPUTATIONS!
   Inputs: xc = matrix of the complete covariates
                (enter 0 if no complete covariates)
           ximp = nimp stacked nonidentical copies of the matrix of incomplete
                  covariates: the copies should be the same at entries that
                  correspond to observed covariate values, but should
                  have potentially different random imputations at entries
                  that correspond to missing covariate values
           y = indicator of case status (1 = case, 0 = control)
           id = vector of id numbers for matched sets
           nc = column vector giving no. observed values for each
                successive column in ximp
           nimp = number of imputations (at least 2)
           bnames = vector of cefficient names
                    (set to 0 if no printed output is desired)
           yname = scalar name of outcome (may be 0 if no output).
 Outputs: bm = mean coefficient estimates
          bmcov = Robins covariance matrix
          devm = mean deviance
          tdf = approximate degrees of freedom for testing each coefficient.
 Global (may be set by user from calling program; default is 0.95):
          _clevel = confidence level (expressed as a proportion)
*/
proc (4) = clrimp(xc,ximp,y,id,nc,nimp,bnames,yname);
declare matrix _clevel = .95;
local t0,nr,ind,k,bcovm,devm,b,bcov,dev,df,bs,bm,bicov,bmcov,tdf;
t0 = date; bnames = 0$+bnames;
if bnames[1];
   "PROC CLRIMP.G: Conditional logistic regression on multiple imputed data.";
endif;
nr = rows(y);
if rows(ximp) ne nr*nimp;
   "INPUT ERROR TO CLRIMP.G: rows(ximp) ne rows(y)*nimp"; end;
   elseif nimp eq 1;
   "INPUT ERROR TO CLRIMP.G: nimp must be 2 or greater."; end;
endif;
nc = sqrt(nr/nc);
if rows(xc) eq nr; nc = ones(cols(xc),1)|nc; else; xc = {}; endif;
/* do nimp regressions: */
   ind = seqa(1,1,nr) - nr;
   k = 0; bs = {}; bcovm = 0; devm = 0; "Regression ";;
   do until k ge nimp; k = k + 1; ind = ind + nr; format 3,0; k;;
      { b,bcov,dev,df } = condlr(xc~ximp[ind,.],y,id,0,0,0);
      bs = bs~b; bcovm = bcovm + bcov; devm = devm + dev;
   endo; print;
/* mean estimates: */ bcovm = bcovm/nimp; bm = meanc(bs'); devm = devm/nimp;
/* between-imputation cov: */
   bicov = ((1 + 1/nimp)/(nimp-1))*moment(bs' - bm',0);
/* with Robins correction: */ bmcov = nc'.*bcovm.*nc + bicov;
/* To set Wald limits, use different t distributions for each bm[j],
   with vector of approximate degrees of freedom: */
   tdf = (nimp - 1)*((diag(bcovm)./(diag(bicov)*(1 +1/nimp)) + 1)^2);
if bnames[1]; local mask4,mask5,fmt4,fmt5,se,p,rad,pll,pul;
   let mask4[1,4] = 0 1 1 1; let mask5[1,5] = 0 1 1 1 1;
   let fmt4[4,3] = "-s" 8 8  "lf" 12 5 "lf," 12 5 "lf" 12 5;
   let fmt5[5,3] = "-s" 8 8  "lf" 12 5 "lf," 12 5 "lf" 17 1 "lf" 9 5;
   se = sqrt(diag(bmcov)); p = 2*cdftc(abs(bm)./se,tdf); format 9,0;
   "Results of regression of ";; $yname;; " on ";; $bnames';; " - ";
   "Estimated betas, se, df, and p: ";
   call printfm(bnames~bm~se~tdf~p,mask5,fmt5);
   rad = invt(_clevel,tdf).*se; pll = exp(bm - rad); pul = exp(bm + rad);
   "Estimated ratios per unit regressors, and ";;
    format 2,0; 100*_clevel;;"% limits:";
   call printfm(bnames~exp(bm)~pll~pul,mask4,fmt4);
   "Mean deviance:";; format 10,3; devm;
   "Total clrimp.g run time: ";; etstr(ethsec(t0,date));
endif;
retp(bm,bmcov,,devm,tdf);
endp;
