/* This procedure fits multiplicative risk ratio models to matched sets using
estimating equations of Greenland (1994), Newton-Raphson fitting
Inputs: y = n-by-2 matrix of set outcomes (exposed pair member in first column)
        x = n-by-p matrix of regressor differences
        rep = n-vector of record repetition counts (1 if none)
        bnames = regressor names
        yname = outcome name (0 for no output)
Locals: b   = vector of estimated coefficients
        bold = previous estimate of beta
        phi	= fitted risk ratios = exp(x*b)
        e   = weighted residuals
        s   = vector of weighted sum of scores
        v   = matrix of negative partial derivatives of s, then inverse
Ouputs: b = coefficients
        v = sandwich covariance estimate for b
        cnv = convergence indicator
Globals: _maxit = maximum number of iterations
         _binit = initial values for b
         _bcriter = convergence criterion for b
*/
// test datasets with one binary regressor, RR = 2
y = {1 1, 1 1, 1 0, 1 0, 1 0, 1 0, 1 0, 1 0, 0 1, 0 1};
x = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
rep = 1; bnames = "helmet"; yname = "death";
{b,v,cnv} = prree(y,x,rep,bnames,yname);
y = {1 1, 1 0, 0 1}; x = {1, 1, 1}; 
let rep = 2 6 2; 
{b,v,cnv} = prree(y,x,rep,bnames,yname);
end;

proc(3) = prree(y,x,rep,bnames,yname);
local iter,bold,phi,e,s,v,cnv;
declare matrix _binit = 0;       // initial values for b
declare matrix _maxit = 30;      // maximum iterations
declare matrix _bcriter = .0005; // convergence criterion for b
format 9,3; bnames = 0$+bnames; yname = 0$+yname;
if yname; print; 
   "PROC RREE.G: Loglinear risk-ratio modeling via estimating equations for matched-pair cohort data";
   if rows(bnames) ne cols(x); "INPUT ERROR: No. names not equal to no. of regressors."; end; endif;
   "No. pairs     : "$+ftocv(sumc(rep.*ones(rows(x),1)),1,0); "No. parameters: "$+ftocv(cols(x),1,0);
print;
endif;
b = _binit; iter=0; 
do until iter gt _maxit; iter = iter+1;
	bold = b;                              // save old b
	if yname; "Iteration "$+ftocv(iter,1,0)$+" -- coefficients:";; b'; endif;
	phi = exp(x*b);	                            // update phi
	e = (y[.,1]-y[.,2].*phi)./(phi+1);          // weighted residuals  
	s = x'(rep.*e);                             // total score
	v = x'((rep.*(y[.,2]+e)./(1 + 1/phi)).*x);  // -score derivative
    trap 1; v = invpd(v);                       // invert v
      if scalerr(v); "SCORE DERIVATIVE NOT POSITIVE DEFINITE IN PREE.G";
         retp(0,0,0); endif; trap 0;
    b = bold + v*s';                       // Newton step
    cnv = prodc(abs(b-bold) .lt _bcriter); // convergence indicator
    if cnv; break; endif;                  // end when all coefficients converge
endo;
if not cnv; "WARNING: no convergence"; endif;
v = v*(x'((rep.*e^2).*x))*v; // sandwich covariance estimate
if yname; call rreport(b,sqrt(diag(v)),bnames,yname,-1,-1,0,1); endif; 
/* if you do not have rreport.g, use:
local se; se = sqrt(diag(v));          // estimated standard errors
"b, estimated SE, Wald p:";; b~se~cdfnt(b./se);
"exp(b) and Wald 95% CL :";; exp(b+(0~(-1.96)~1.96)*se); */
retp(b,v,cnv);
endp;
