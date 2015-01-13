/* This procedure does unconditional logistic regression
   on multiple-imputed data using Rubin's method.
   A stack of multiple-imputed covariate matrices suitable for input to this
   procedure (denoted by ximp below) may be generated from MNIMP.G.
   Inputs: xc = matrix of the complete covariates
                (enter 0 if no complete covariates)
           ximp = nimp stacked nonidentical copies of the matrix of incomplete
                  covariates: the copies should be the same at entries that
                  correspond to observed covariate values, but should
                  have potentially different random imputations at entries
                  that correspond to missing covariate values
           y = indicator of case status (1 = case, 0 = noncase)
           nc = column vector giving no. observed values for each
                successive column in ximp
           nimp = number of imputations (at least 2)
           const = 1 if constant (intercept) term to be added, 0 if not
           bnames = vector of coefficient names
                  (set to 0 if no printed output is desired)
           yname = scalar name of outcome (may be 0 if no output).
 Outputs: bm = mean coefficient estimates
          bmcov = Robins covariance matrix
          devm = mean deviance
          tdf = approximate degrees of freedom for testing each coefficient.
 Global (may be set by user from calling program; default is 0.95):
          _clevel = confidence level (expressed as a proportion)
*/
proc (4) = lrimp(xc,ximp,y,nc,nimp,const,bnames,yname);
declare matrix _clevel = .95;
local t0,nr,np,ind,k,bs,bcovm,devm,b,bcov,bcovr,dev,rss,df,bm,bicov,bmcov,tdf;
t0 = date; bnames = 0$+bnames;
"PROC LRIMP.G: Logistic regression on multiple imputed data.";
nr = rows(y);
if rows(ximp) ne nr*nimp;
   "INPUT ERROR TO LRIMP.G: rows(ximp) ne rows(y)*nimp"; end;
   elseif nimp eq 1;
   "INPUT ERROR TO LRIMP.G: nimp must be 2 or greater."; end;
endif;
if const; xc = ones(nr,1)~xc; bnames = "constant"|bnames; endif;
nc = sqrt(nr/nc);
if rows(xc) eq nr; nc = ones(cols(xc),1)|nc; else; xc = {}; endif;
/* do nimp regressions: */
   ind = seqa(1,1,nr) - nr;
   k = 0; bs = {}; bcovm = 0; devm = 0; "Imputation ";;
   do until k ge nimp; k = k + 1; ind = ind + nr; format 3,0; k;;
      { b,bcov,bcovr,dev,rss,df } = logreg(xc~ximp[ind,.],y,1,0,0,0,0,0);
      bs = bs~b; bcovm = bcovm + bcov; devm = devm + dev;
   endo; print;
/* mean estimates: */ bcovm = bcovm/nimp; bm = meanc(bs'); devm = devm/nimp;
/* between-imputation cov: */
   bicov = ((1 + 1/nimp)/(nimp-1))*moment(bs' - bm',0);
/* with Robins correction: */ bmcov = nc'.*bcovm.*nc + bicov;
/* To set Wald limits, use different t distributions for each bm[j],
   with vector of approximate degrees of freedom: */
   tdf = (nimp - 1)*((diag(bcovm)./(diag(bicov)*(1 +1/nimp)) + 1)^2);
if bnames[1+const]; local mask4,mask5,fmt4,fmt5,se,p,rad,pll,pul,pr;
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
   pr = bnames~exp(bm)~pll~pul;
   call printfm(pr[1+const:rows(bm),.],mask4,fmt4);
   "Mean deviance:";; format 10,3; devm;
   "Total lrimp.g run time: ";; etstr(ethsec(t0,date));
endif;
retp(bm,bmcov,devm,tdf);
endp;
