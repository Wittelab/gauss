/* Test program:  let x = -1 0 1; let y = -3 1 2; let w = 1 2 1;
   let rep = 20 10 20; rep = 1; const = 1; x = x~(x^2);
   { b,bcov,rss,rsq,df } = wls(x,y,w,rep,const,"x"|"x^2","y");
   xf = 2*sortc(rndu(50,1),1) - 1; xf = xf~(xf^2); _pdate = 0;
   r = 1.96; call ypred(xf,b,bcov,const,r,0,2,"x","y");
   end; */
/* This proc computes the fitted mean and confidence intervals in a GLM
   at the regressor points (rows) in x from the fitted coefficients b with
   covariance matrix c.
   Inputs: xf = regressor points at which to compute fitted mean
           b = coefficient estimates
           bcov = estimated covariance matrix for b
           const = 1 if constant to be added to x, 0 if not
           r = multiplier for confidence percentile, e.g., 1.960 for a 95%
               normal interval, 2.228 for a 95% t interval on 10 d.f.
           link = code for link function:
                  Enter 0 for identity (linear model),
                  1 for log-linear models, 2 for logistic
                  (anything other than 1 or 2 will use identity).
           g = index of column of xf against which to plot fitted means
               (0 for no graph)
           xname = name of xf[.,g] (may be 0 if no graph)
           yname = name of regressand y (may be 0 if no graph)
   Output: yhat = 3-column matrix with the fitted means,
                  lower limits, and upper limits for each row of x.
*/
proc ypred(xf,b,bcov,const,r,link,g,xname,yname);
local yhat;
if const; xf = ones(rows(xf),1)~xf; endif;
yhat = xf*b + (0~(-1)~1).*(r*sqrt(sumr((xf*bcov).*xf)));
if link eq 1; yhat = exp(yhat); elseif link eq 2; yhat = expit(yhat); endif;
if g; library pgraph;  /* SET GRAPHICS PARAMETERS -- */ graphset;
    xlabel(xname); ylabel(yname);
    /* line types: */ let _pltype =  6 4 4;
    /* line widths: */ let _plwidth =  1 3 3;
    call xy(xf[.,g],yhat);
endif;
retp(yhat);
endp;

