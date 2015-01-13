/* Examples for LRMSM.G:  output file = lrmsm.out reset;
let x = 1 1 0 0 1 1 0 0; z0 = 0; let z = 1 1 1 1 0 0 0 0;
let y = 1 0 1 0 1 0 1 0; t = 0; id = seqa(1,1,8);
aname = "x"; let xnames = "x" "z"; yname = "disease"; _collaps = 1;
/* defaults: */ _current = 1; _alag = -1; _ylag = -1; _rlag = -1; _mis = 0;
let rep = 80 20 60 40 40 60 30 120;
call logreg(x,y,1,rep,1,0,xnames[1],yname);
ian = 0; iad = 2; iy = 1;
call lrmsm(x,x~z,y,0,t,id,ian,iad,iy,0,0,rep,1,0,aname,xnames,yname);
let rep = 80 20 60 40 40 60 20 80;
call logreg(x,y,1,rep,1,0,xnames[1],yname);
call lrmsm(x,x~z,y,0,t,id,ian,iad,iy,0,0,rep,1,0,aname,xnames,yname);
end; */
/* This procedure fits marginal structural logistic models to a polytomous and
   possibly time-dependent treatment x and a binary outcome y.
   This involves three regressions, plus two more if data are missing.
   NOTE 1: The input data must be arranged with a separate record for each
           subject at each time observed for that subject. Thus, a subject
           observed at k times would have k records.
   NOTE 2: The records for each subject must be contiguous,
           and sorted from earliest to latest time.
   NOTE 3: For valid effect estimation, all time-dependent covariates other
           than treatment history must be entered as marginalizing covariates;
           equivalently, all conditioning covariates except treatment history
           must be fixed (constant within subjects).
   NOTE 4: The repetition count must apply to subjects (clusters), not records,
           and so should be constant within subjects.
   Inputs: a = vector of current-treatment indicators
           x = matrix containing all covariates for all regressions
           y = vector of outcome indicators
           r = vector of indicators of missing current data (0 if none)
           t = time vector (times must be nonnegative)
           id = subject ID no.
           ian = indices of columns in x containing treatment-history and fixed
                 covariates for current-treatment regression (0 if none)
           iad = indices of columns in x containing all covariates for
                 current-treatment regression
           iy = indices of columns in x containing covariates for outcome
                regression
           irn = indices of columns in x containing treatment-history and fixed
                 covariates for missing-data regression (0 if none)
           ird = indices of columns in x containing all covariates for
                 missing-data regression (0 if none)
           rep = repetition count (0 if none)
                 --> MUST BE CONSTANT WITHIN SUBJECTS (cluster repetition)
           const = 1 if program to add constant, 0 if not
           offset = offset vector (0 if none)
           anames = name of a column (0 for no printed output)
           xnames = names of x columns (0 for no printed output)
           ynames = name of outcome (0 for no printed output)
   Outputs: b = coefficient estimates (with intercept first if const = 1)
            bcov = robust ("sandwich") covariance matrix for b,
                   corrected for estimation of weights from treatment and
                   missing-data regression
            bcovr = uncorrected robust covariance matrix for b
            dev = residual weighted deviance for model
            rss = residual weighted sum of squares
            df = residual degrees of freedom
Globals (may be set by user from calling program; otherwise defaults are used):
      _current = 1 to use current exposure as a regressor for y,
                 0 if not (default is 1)
      _alag = lag time for a (records with times le _alag not used in
              treatment regression -- default is -1)
      _ylag = lag time for y (records with times le _ylag not used in
              outcome regression -- default is -1)
      _rlag = lag time for r (records with times le _rlag not used in
              missing-data regression -- default is -1)
      _clevel = confidence level expressed as a proportion (default is .95)
      _fit = 1 for iteration reporting and tests of fit, 0 if not (default is
             1, but no tests provided if mean of case or noncase count is <3)
      _mis = 0 if no missing values, 1 if missing data
      _mcode = missing-value scalar or column vector with element for each
               column of x~z~z0~y (default -99999; ignored if _mis = 0)
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _maxit = maximum number of iterations (default is 30)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
*/
proc (6) =
     lrmsm(a,x,y,r,t,id,ian,iad,iy,irn,ird,rep,const,offset,aname,xnames,yname);
declare matrix _current = 1; declare matrix _alag = -1;
declare matrix _xlag = -1; declare matrix _ylag = -1;
declare matrix _mis = 0; declare matrix _amcode = -99999;
declare matrix _xmcode = -99999; declare matrix _ymcode = -99999;
declare matrix _binit = 0; declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .1;
declare matrix _fit = 1; declare matrix _clevel = .95;
local t0,nr,ra,rx,mcc,c,m,b,bcov,dev,rss,df,p,nmu,dmu,rscore,rw,
      ascore,aw,u,bscore,bcovr;
t0 = date;
aname = 0$+aname; xnames = 0$+xnames; yname = 0$+yname;
nr = rows(x); c = ones(nr,1);
if _mis; /* indicate nonmissing data: */
   ra = a .ne _amcode'; rx = x .ne _xmcode'; r = (y .ne _ymcode).*r;
   else; ra = c; rx = ones(nr,cols(x));
endif;
if rows(rep) eq 1; rep = c; endif;
if xnames[1]; print; format /rds 8,2;
   if rows(aname) ne cols(a) or rows(xnames) ne cols(x);
      "INPUT ERROR: No. names supplied not equal to no. of covariates."; end;
   endif;
   "PROC LRMSM.G: Marginal structural logistic regression of "$+yname$+" on";
   $xnames[iy]';
   if sumc(iad); " Treatment predictors: "; $xnames[iad]';
      else; " No treatment predictors."; endif;
   if sumc(ird); " Missingness predictors: "; $xnames[ird]';
      else; " No missingness predictors."; endif;
endif;
if rank(x[.,iy]'x[.,iy]) lt cols(x[.,iy]);
   "DESIGN MATRIX RANK DEFICIENT IN LRMSM.G"; retp(0,0,0,0,0,0); endif;
if sumc(ird);
   if xnames[1]; print;
      " *******> REGRESSING NONMISSINGNESS ON COVARIATES <*******"; endif;
   m = (t .gt _rlag).*prodr(rx[.,ird]);
   if sumc(irn);
     { b,bcov,bcovr,dev,rss,df } =
       logreg(c~x[.,irn],r,1,m.*rep,0,0,"const"|xnames[irn],"nonmiss");
     p = exp((c~x[.,irn])*b);
     else; p = m'r; p = p./(sumc(m) - p);
   endif;
   nmu = (p^r)./(1 + sumr(p));
   { b,bcov,bcovr,dev,rss,df } =
     logreg(c~x[.,ird],r,1,m.*rep,0,0,"const"|xnames[ird],"nonmiss");
   p = exp((c~x[.,ird])*b);
   dmu = (p^r)./(1 + sumr(p));
   rscore = matcol(id~(x[.,ird].*(r-dmu)),1);
   rscore = rscore[.,2:cols(rscore)];
   rw = nmu./dmu; rw = m.*rw + (1-m);
   /* Multiply weights over time periods: */
      u = cprodg(id~rw,1); rw = u[.,2:cols(u)];
   else; rscore = {}; r = c; rw = r;
endif;
if sumc(iad);
   if xnames[1]; print;
      " *******> REGRESSING TREATMENT ON COVARIATES <*******"; endif;
   m = (t .gt _alag).*prodr(ra~rx[.,iad]);
   if sumc(ian);
      { b,bcov,dev,df } =
        lrpoly(c~x[.,ian],a,1,m.*rep,0,0,"const"|xnames[ian],aname);
      p = exp((c~x[.,ian])*b);
      else; p = m'a; p = p./(sumc(m) - p);
   endif;
   nmu = sumr((a.*p)~(1-sumr(a)))./(1 + sumr(p));
   { b,bcov,dev,df } =
     lrpoly(c~x[.,iad],a,1,m.*rep,0,0,"const"|xnames[iad],aname);
   p = exp((c~x[.,iad])*b);
   dmu = sumr((a.*p)~(1-sumr(a)))./(1 + sumr(p));
   ascore = matcol(id~(x[.,iad].*(a-dmu)),1);
   ascore = ascore[.,2:cols(ascore)];
   aw = nmu./dmu; aw = m.*aw + (1-m);
   /* Multiply weights over time periods: */
      u = cprodg(id~aw,1); aw = u[.,2:cols(u)];
   if not _current; aw = aw.*dmu./nmu; endif;
   else; rscore = {}; aw = c;
endif;
/* Multiply weights together: */
   u = rw.*aw.*(t .gt _ylag).*prodr(r~rx[.,iy]);

print; " *******> FINAL STEP: FIT MSM TO OUTCOME <*******"; print;

if xnames[1]; mcc = meancw(rep.*(y~(1-y)),u .gt 0);
   "No. of cases and noncases used  : ";; format 7,0; ((u .gt 0)'rep)*mcc';
   "Minimum and maximum used weights: ";;
    format 5,3; minc(selif(u,u .gt 0));; maxc(u);
   "Weights will be geometrically recentered for outcome regression.";
endif;
u = u./exp(meanc(ln(selif(u,u .gt 0))));
x = x[.,iy];
{ b,bcov,bcovr,dev,rss,df } = logreg(x,y,1,u.*rep,const,0,0,0);
if const eq 1; x = c~x; endif;
bscore = matcol(id~rep~(x.*(u.*(y-expit(x*b)))),1|2);
rep = bscore[.,2]; bscore = bscore[.,3:cols(bscore)];
bcovr = bcov*((bscore.*rep)'bscore)*bcov;
if sumc(ird|iad); /* corrected robust covariance: */
   local bra,bracov,bar,r2,bradf;
   { bra,bracov,bar,r2,bradf } = wls(rscore~ascore,bscore,1,rep,0,0,0);
   bscore = bscore - (rscore~ascore)*bra;
   bcov = bcov*((bscore.*rep)'bscore)*bcov;
endif;
if xnames[1]; print;
   "Using uncorrected specification-robust covariance matrix:";
   rreport(b,sqrt(diag(bcovr)),xnames[iy],yname,-1,
           rss*(-1)^(1-_fit),df,1+_fit*(minc(mcc) ge 3));
   if sumc(ird|iad); print;
      "Using corrected specification-robust covariance matrix:";
      rreport(b,sqrt(diag(bcov)),xnames[iy],yname,-1,-1,df,2);
   endif;
   if _fit;
      "The means of the case & noncase counts are ";; format 6,1; mcc';
      if minc(mcc) lt 5;
         "The minimum of the mean case & mean noncase counts is less than 5.";
         if minc(mcc) ge 3;
           "Hence the above test of fit may be invalid.";
         endif;
      endif;
      print; "Total LRMSM.G run time: ";; etstr(ethsec(t0,date));
   endif;
   print; format 10,3;
endif;
retp(b,bcov,bcovr,dev,rss,df);
endp;
