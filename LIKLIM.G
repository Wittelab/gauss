/* This proc searches for LR limits.
 Inputs: x = design matrix,
         y = vector of case counts,
         n = vector of totals,
         rep = vector of repetion counts (0 if none),
         const = 1 if constant must be added to x, O otherwise
                 (limits won't be computed for constant)
         offsetin = vector of offsets (0 if no offsets),
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired)
         yname = outcome name
         indl = k-vector of indices of bhat for which limits will be found
                (all if 0)
  Outputs: limits = k by 2 matrix of lower and upper limits for b[indl].
Globals (may be set by user from calling program; otherwise defaults are used):
      _clevel = confidence level expressed as a proportion (default is .95)
      _mis = 0 if no missing values, 1 if complete-case analysis
      _mcode = missing-value scalar or column vector with element
               for each column of x (default is -99999; ignored if _mis = 0)
      _bcriter = convergence criterion for b (default is .001) (for lrquickl)
      _criter = convergence criterion for limits
      _report = 1 to report results at each iteration (default is 0)
*/
proc liklim(x,y,n,rep,const,offsetin,bnames,yname,indl);
/* confidence level: */ declare _clevel = .95;
/* zero initial values: */ declare matrix _binit = 0;
/* convergence criterion: */ declare matrix _criter = .001;
/* max. iterations: */ declare matrix _maxit = 30;
declare matrix _report = 0;
local t0,np,bhat,cnvl,se,z,limits,j,ind,irs,i,cnv,rad,offset,b,bcov,
      llmax,llmod,rsst,df,radold;
t0 = date; bnames = 0$+bnames;
if bnames[1]; print;
   " Proc liklim.g: likelihood limits for logistic regression."; endif;
n = n.*ones(rows(x),1);
if indl[1] eq 0; indl = seqa(1,1,cols(x)); endif;
if sumc(rep) ne 0; y = y.*rep; n = n.*rep; endif;
if const; /* Include constant if requested: */
   x = ones(rows(x),1)~x;
   indl = 1 + indl; bnames = "constant"|bnames;
endif;
np = cols(x);
if rows(_binit) eq 1; bhat = zeros(np,1); else; bhat = _binit; endif;
{ bhat,bcov,llmax,cnvl } = lrquickl(x,y,n,offsetin,bhat);
se = sqrt(diag(bcov));
/* Desired Z value: */ z = cinvnorm((1-_clevel)/2);
/* Storage for limits: */ limits = zeros(rows(indl),2);
/* Loop over indl: */ j = 0;
do until j ge rows(indl); j = j + 1;
   /* Index of variable to be offset (for which limits are
      to be tested): */ ind = indl[j];
   /* Define vector of regressor indices (excluding ind): */
      if np eq 1; irs = 1;
         elseif ind eq 1 /* ind is first variable */; irs = seqa(2,1,np-1);
         elseif ind eq np /* ind is last variable */; irs = seqa(1,1,np-1);
         else; irs = seqa(1,1,ind-1)|seqa(ind+1,1,np-ind);
      endif;
   /* Initialize CI radius at Wald radius: */ rad = z*se[ind];
   if bnames[1] /* print headers for this variable: */;
      format /rdn 10,5; print;
      " ML OR estimate & Wald limits for "$+bnames[ind]$+":";;
        exp(bhat[ind]+(0~(-rad)~rad)); format 3,0;
      100*(1-cdfnt(z));; format 10,5;
      "% likelihood limits for "$+bnames[ind]$+":";
   endif;
   /* Iterate to lower limit using offset to find where deviance test yields
      desired z value: */
   i = 0; cnv = 0; b = bhat[irs];
   do until (i gt _maxit) or cnv; i = i + 1;
      limits[j,1] = real(bhat[ind] - rad);
      offset = limits[j,1]*x[.,ind];
      if np eq 1;
         llmod = y'(offset)+n'ln(n.*expit(-offset)+(n.*expit(-offset) .eq 0));
         else; {b,bcov,llmod,cnvl} = lrquickl(x[.,irs],y,n,offset+offsetin,b);
               if cnvl le 0; break; endif; b = real(b);
      endif;
      radold = rad;
      rad = rad*(z/sqrt(2*abs(llmax-llmod)))^sign(llmax-llmod);
      /* Test for convergence: */ cnv = abs(rad - radold) lt _criter;
      if _report or bnames[1]*(i eq 1);
        format /rdn 3,0; "     At iteration ";; i;; format 10,5;
        ", lower limit & p = ";; exp(limits[j,1]);; cdfchic(2*(llmax-llmod),1);
      endif;
   endo;
   if bnames[1] and cnvl ge 0;
      "   After "$+ftocv(i,1,0)$+" iterations, lower OR limit & deviance p:";;
      exp(limits[j,1]);; cdfchic(2*(llmax-llmod),1);
      else; "Failed to find lower limit";
   endif;
   rad = z*se[ind];
   /* Iterate to upper limit: */
   i = 0; cnv = 0; b = bhat[irs];
   do until (i gt _maxit) or cnv; i = i + 1;
      limits[j,2] = real(bhat[ind] + rad);
      offset = limits[j,2]*x[.,ind];
      if np eq 1;
         llmod = y'(offset)+n'ln(n.*expit(-offset)+(n.*expit(-offset) .eq 0));
         else; {b,bcov,llmod,cnvl} = lrquickl(x[.,irs],y,n,offset+offsetin,b);
               if cnvl le 0; break; endif; b = real(b);
      endif;
      radold = rad;
      rad = rad*(z/sqrt(2*abs(llmax-llmod)))^sign(llmax-llmod);
      cnv = abs(rad - radold) lt _criter;
      if _report or bnames[1]*(i eq 1);
        format /rdn 3,0; "     At iteration ";; i;; format 10,5;
        ", upper limit & p = ";; exp(limits[j,2]);; cdfchic(2*(llmax-llmod),1);
      endif;
   endo;
   if bnames[1] and cnvl gt 0;
      "   After "$+ftocv(i,1,0)$+" iterations, upper OR limit & deviance p:";;
      exp(limits[j,2]);; cdfchic(2*(llmax-llmod),1);
      else; "Failed to find upper limit";
   endif;
endo;
if bnames[1];"Total run time: ";;etstr(ethsec(t0,date));format /rds 9,3;endif;
retp(limits);
endp;

