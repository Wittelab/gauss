/* This procedure does ordinal ML logistic regression using the
   continuation-ratio model, via nested ordinary logistic regression.
 Inputs: x = design matrix,
         y = matrix of outcome counts (do not include reference level),
         n = vector of totals,
         rep = vector of repetion counts (0 if none),
         offset = vector of offsets (0 if no offsets);
                  matrix of offsets may be used if cols(offset) = cols(y),
                  in which case a separate offset willbe used for each outcome.
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
         ynames = vector of outcome names (may be 0 if no printed output).
Outputs: b = intercept estimates (one per index category)
             and slope coefficient estimates,
         bcov = inverse-information covariance estimate for b,
         bcovr = robust (sandwich) covariance estimate for b,
         dev = residual deviance for model,
         rss = residual sum of squares for model,
         df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _maxit = maximum number of iterations (default is 30)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
      _specrob = 1 if you want an additional version of output
                 in which the Huber specification-robust "sandwich"
                 covariance matrix estimate (bcovr) is used instead of
                 the usual inverse-information estimate (default is 0)
*/
proc (6) = lrordcr(x,y,n,rep,offset,bnames,yname);
/* _binit and _specrob are used only by logreg.g, which is called at end. */
declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .1;
local nd,c,b,bcov,bcovr,dev,rss,df;
bnames = 0$+bnames;
if sumc(rep); y = y.*rep; n = n.*rep; endif;
/* No. index outcome categories: */ nd = cols(y);
if bnames[1]; print; format /rdn 3,0;
   "Proc lrordcr.g: Continuation-ratio model with";; nd;; " index outcomes";
   "     labelled y = ";; nd;; " to y = 1, and ";; cols(x);;" regressors.";
   if rows(bnames) ne cols(x);
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
endif;
/* Rearrange data for logistic regression by 
   1) making a matrix of outcome-category indicators; 
   2) stacking nd copies of x on top of each other;
   3) summing rows of y and n-sumr(y) to make nested denominators; and
   4) summing rows of y to make nested numerators: */
/* 1. Make a separate intercept for each of the nd index levels of y: */
   c = eye(nd).*.ones(rows(x),1);
/* 2. Stack nd copies of x on top of each other: */ x = ones(nd,1).*.x;
if rows(offset) gt 1; /* stack nd copies of offset on top of each other: */
   offset = ones(nd,1).*.offset;
endif;
/* 3. Make new, nested n in two steps - 
   Step A:
   rev(y') is the matrix in which row j is the counts at index level nd+1-j.
   The row of counts at the reference level are (n - sumr(y))'.
   cumsumc creates a matrix of partial sums down the rows.
   n thus becomes the matrix in which column j is the total number of subjects
   at or above level nd+1-j of y: */ 
   n = cumsumc(rev(y')|(n-sumr(y))')';
/* B: These nested denominators are next "rolled out" without the first column,
   to become the vector of denominators for the nested logistic regression: */
   n = vec(n[.,2:(nd+1)]);
/* 4. Make new, nested y -
   cumsumc(rev(y')) is the matrix in which column j gives the counts above
   level nd+1-j of y; this matrix is rolled out to become the vector of
   "cases" for the nested logistic regression: */
   y = vec(cumsumc(rev(y'))');
if bnames[1]; /* create names for the nd constants: */
   bnames = (0$+"y = "$+ftocv(seqa(nd,-1,nd),2,0))|bnames;
endif;
/* Run the nested logistic regression: */
{ b,bcov,bcovr,dev,rss,df } = logreg(c~x,y,n,0,0,offset,bnames,yname);
retp(b,bcov,bcovr,dev,rss,df);
endp;
