/* This program does conditional logistic regression for m:(n-m) matched sets
   with variable matching ratio.
   Inputs: x = regressor matrix
           y = indicator of case status (1 = case, 0 = noncase)
           id = vector of ID numbers of matched sets to which subjects belong
           reps = repetion count for sets (set to 0 if none) by SORTED id no.
                (reps[i] must refer to the count for the i'th ID no.)
           Note that rows(id) = rows(x) > rows(reps).
           offset = vector of offsets (set to zero if no offsets)
           mprin = 1 if you want description of matched sets, 0 otherwise
           bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
           yname = scalar name of outcome (may be 0 if no output).
   Outputs: b = coefficient estimates,
            bcov = inverse-information covariance matrix,
            dev = -2*conditional loglikelihood for model,
            np = model degrees of freedom (no. parameters)
Globals (may be set by user from calling program; otherwise defaults are used):
      _mis = 0 if no missing values, 1 if complete-case analysis
      _mcode = missing-value scalar or column vector with element
               for each column of x (default -99999; ignored if _mis = 0)
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _maxit = maximum number of iterations (default is 30)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
      _fit = 1 to print fitting and global test results (default is 0)
*/
proc (4) = condlrm(x,y,id,reps,offset,mprin,bnames,yname);
declare matrix _mis = 0; declare matrix _mcode = -99999;
declare matrix _binit = 0; declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .1;
declare matrix _fit = 1;
local t0,np,ind,uid,del,n,sn,m,i,xs,ns,ids,ofs,xi,yi;
t0 = date; bnames = 0$+bnames; yname = 0$+yname;
/* No. parameters: */ np = cols(x);
if rows(offset) eq 1; offset = offset*ones(rows(x),1); endif;
/* vector of sorted indices: */ ind = sortind(id);
/* sorted covariates matrix: */ x = x[ind,.];
/* sorted outcome vector: */ y = y[ind];
/* sorted id nos.: */ id = id[ind];
/* sorted offsets: */ offset = offset[ind];
/* vector of unique id nos.: */ uid = unique(id,1);
if rows(reps) eq 1; /* set rep to 1: */ reps = ones(rows(uid),1); endif;
if _mis eq 1; /* delete records with missing regressor values: */
   del = sumr((x~y~id~offset) .eq _mcode') .gt 0;
   if sumc(del); x = delif(x,del); y = delif(y,del); id = delif(id,del);
      offset = delif(offset,del); del = unique(id,1);
      if rows(del) lt rows(uid);
         reps = selif(reps,sumr(uid .eq del')); uid = del; endif;
   endif;
endif;
/* For each matched set, compute the size & the no. of cases: */
   /* total no. in each matched set: */ n = counts(id,uid);
   /* no. of cases in each matched set: */ m = countwts(id,uid,y);
if bnames[1]; print;
   "PROC CONDLRM.G: conditional logistic regression."; print;
   "No. of cases and controls: ";; format /rds 5,0; sumc(reps.*(m~(n-m)))';
   if mprin; /* print description of the sets: */
      "Matched set ID no., no. cases, no. noncases, & rep. count:";
      uid~m~(n-m)~reps; print;
   endif;
   if rows(bnames) ne cols(x);
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
endif;
/* Eliminate sets with no case or no control: */ del = m.*(m-n) .eq 0;
   if sumc(del); n = delif(n,del); m = delif(m,del); uid = delif(uid,del);
      reps = delif(reps,del); del = 1-sumr(id .eq uid');
      x = delif(x,del); y = delif(y,del);
      id = delif(id,del); offset = delif(offset,del);
   endif;
if (_mis or sumc(del)) and bnames[1];
   "No. of cases and controls after eliminating sets";
   "with no complete case or no complete control record:";;
   sumc(reps.*(m~(n-m)))';
endif;
if rank(x) lt np; "DESIGN MATRIX RANK DEFICIENT IN CONDLRM.G";
   retp(0,0,0,0); endif;
/* Construct a new design matrix xs, consisting of covariate sums of all
   possible subsets of size m from the n subjects in each matched set,
   minus case covariate sums, for subsets containing at least one control: */
   /* positions in x where each matched set starts: */
      sn = cumsumc(n) - n + 1;
   /* loop through the sets: */
      i = 0; ns = {}; ids = {}; xs = {}; ofs = {};
   do until i ge rows(n); i = i + 1;
      /* positions of set i members in x: */ ind = seqa(sn[i],1,n[i]);
      /* Covariate and outcome sums for all subsets of size m[i]: */
         xi = sumcomb(y[ind]~offset[ind]~x[ind,.],m[i]);
      /* Remove subsets with cases only: */ xi = selif(xi,xi[.,1] .lt m[i]);
      /* No. of subsets left: */ ns = ns|rows(xi);
      /* Matched-set id numbers for rows of xs: */
         ids = ids|(i*ones(rows(xi),1));
      /* Subtract off the case covariate totals, then paste xi to xs: */
         xs = xs|(xi[.,3:np+2] - y[ind]'x[ind,.]);
         ofs = ofs|(xi[.,2] - y[ind]'offset[ind]);
   endo;
/* Begin ML logistic fitting of constant outcome to xs, with no intercept: */
   local b,bold,bcov,dev,dev0,devold,iter,cnv,score,inf,ss;
   /* Note that saturated loglikelhood = 0 */
   /* Initialize beta, deviance, counter, convergence indicator: */
      if rows(_binit) eq np; b = _binit;
         elseif sumc(reps) eq rows(reps);
         /* start with unmatched WLS estimates: */ local p,w,u;
         p = (y + sumc(y)/rows(y))/2; w = p.*(1-p);
         u = ones(rows(x),1);
         trap 1; bcov = invpd(moment(sqrt(w).*(u~x),0));
         if scalerr(bcov); "SINGULARITY INITIALIZING CONDLRM.G";
            retp(0,0,0,0); endif; trap 0;
         b = bcov*((u~x)'(w.*(logit(p)-offset))); b = b[2:rows(b)];
         else; b = zeros(np,1);
      endif;
      dev = 0; iter = 0; cnv = 0;
   do until cnv or (iter ge _maxit); iter = iter + 1;
      /* save old values: */ devold = dev; bold = b;
      /* compute conditional deviance, score, & information at b: */
         { dev,score,inf } = mscore(xs,ofs,b,ns,reps,ids);
      /* covariance of b: */ trap 1; bcov = invpd(inf);
         if scalerr(bcov); "SINGULARITY IN CONDLRM.G";
            retp(0,0,0,0); endif; trap 0;
      /* update b: */ b = b + bcov*score;
      if _fit and bnames[1];
         "Deviance at iteration ";; format 4,0; iter;;":";; format 11,2; dev;
      endif;
      cnv = prodc(abs(b-bold) .lt _bcriter) and (abs(devold-dev) lt _dcriter);
   endo;
  if iter gt _maxit;
     "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
  endif;
if bnames[1]; /* print results: */
   rreport(b,sqrt(diag(bcov)),bnames,yname,-1,-1,0,1);
   if _fit;
     "Model degrees of freedom:";; format 5,0; np;
     { dev0,score,inf } = mscore(xs,ofs,zeros(np,1),ns,reps,ids);
     "Deviance statistic and p for all coefficients:";;
     format 11,3; dev0-dev;; format 11,5; cdfchic(dev0-dev,np);
     "Score statistic and p for all coefficients   :";;
     ss = score'invpd(inf)*score;
     format 11,3; ss;; format 11,5; cdfchic(ss,np);
     "Total run time: ";; etstr(ethsec(t0,date));
   endif;
endif;
retp(b,bcov,dev,np);
endp;
