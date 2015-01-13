/* This program does conditional logistic regression for 1:m matched sets
   with variable matching ratios.
   WARNINGS:1) EACH MATCHED SET IS ASSUMED TO HAVE AT LEAST ONE CASE AND
            ONE NONCASE.
            2) NO MISSING VALUES OR OFFSETS ARE ALLOWED.
            3) THE INPUT VECTORS AND MATRICES ARE ASSUMED TO BE SORTED ON id,
            THE VECTOR OF MATCHED-SET ID NUMBERS, so that all subjects
            in a matched set should be contiguous in the analysis data.
   CONDLRM.G DOES NOT HAVE THESE RESTRICTIONS.
   To sort your input data on id, you may use this code:
         /* vector of sorted indices: */ ind = sortind(id);
         /* sorted covariates matrix: */ x = x[ind,.];
         /* sorted outcome vector: */ y = y[ind];
         /* sorted id nos.: */ id = id[ind];
         /* sorted rep counts (if present): */ reps = reps[ind];
   Inputs: x = regressor matrix
           y = indicator of case status (1 = case, 0 = control)
           id = vector of id numbers of matched sets from smallest to largest
           reps = repetion count for sets (set to 0 if none) by sorted id no.
                  (reps[i] must refer to the count for the i'th set)
           Note that rows(id) = rows(x) > rows(reps).
           bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
           yname = scalar name of outcome (may be 0 if no output).
   Outputs: b = coefficient estimates,
            bcov = inverse-information covariance matrix,
            dev = -2*conditional loglikelihood for model,
            np = model degrees of freedom (no. parameters)
Globals (may be set by user from calling program; otherwise defaults are used):
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _maxit = maximum number of iterations (default is 30)
               -- if 0, there will be no updating of _binit
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
*/
proc (4) = condlr(x,y,id,reps,bnames,yname);
declare matrix _binit = 0; declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .1;
local t0,uid,m,xs,b,bold,bcov,dev,dev0,nsets,
      devold,iter,cnv,ids,score,inf,np,ss;
t0 = date; bnames = 0$+bnames; yname = 0$+yname;
uid = unique(id,1);
if rows(uid) ne sumc(y);
   "INPUT ERROR IN CONDLR.G:"; format /rds 3,0;
   " Number of unique ID numbers, ";; rows(uid);
   " not equal to number of unique case records, ";;
   sumc(y); format 9,3; end;
endif;
if rows(uid) gt sumc(1-y);
   "INPUT ERROR IN CONDLR.G:"; format /rds 3,0;
   " Number of unique ID numbers, ";; rows(uid);
   " greater than number of unique noncase records, ";;
   sumc(1-y); format 9,3; end;
endif;
/* No. parameters: */ np = cols(x);
if rank(x) lt np; "DESIGN MATRIX RANK DEFICIENT IN CONDLR.G";
   retp(0,0,0,0); endif;
/* sorted control counts: */ m = counts(delif(id,y),uid);
if rows(reps) eq 1; reps = ones(rows(m),1);
   elseif rows(reps) ne rows(m);
   "INPUT ERROR: Number of rep counts must equal number of matched sets."; end;
endif;
if bnames[1]; print;
   "PROC CONDLR.G: 1:m conditional logistic regression.";
   "Number of cases and controls:";; format /rds 5,0; sumc(reps);; reps'm;
   if rows(bnames) ne cols(x);
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
endif;
/* Control minus matched-case data: */
   xs = delif(x,y) - mexpand(selif(x,y),m);
/* Indices for mscore: */
   ids = mexpand(seqa(1,1,rows(m)),m);
/* Note that saturated loglikelhood = 0 */
/* Initialize beta, deviance, counter, convergence indicator: */
   if rows(_binit) eq np; b = _binit;
      elseif sumc(reps) eq rows(reps);
      /* start with unmatched WLS estimates: */ local p,w,u;
      p = (y + sumc(y)/rows(y))/2; w = p.*(1-p);
      u = ones(rows(x),1);
      trap 1; bcov = invpd(moment(sqrt(w).*(u~x),0));
      if scalerr(bcov); "SINGULARITY INITIALIZING CONDLR.G";
         retp(0,0,0,0); endif; trap 0;
       b = bcov*((u~x)'(w.*logit(p))); b = b[2:rows(b)];
       else; b = zeros(np,1);
   endif;
   dev = 0; iter = 0; cnv = 0;
do until cnv or (iter ge _maxit and _maxit ne 0); iter = iter + 1;
   /* save old values: */ devold = dev; bold = b;
   /* compute conditional deviance,score,information at b: */
   { dev,score,inf } = mscore(xs,0,b,m,reps,ids);
   /* covariance of b: */ trap 1; bcov = invpd(inf);
      if scalerr(bcov); "SINGULARITY IN CONDLR.G";
         retp(0,0,0,0); endif; trap 0;
   if _maxit eq 0; break; endif;
   /* update b: */ b = b + bcov*score;
   if bnames[1]; "Deviance at iteration ";;
      format /rds 4,0; iter;; ":";; format 11,2; dev; endif;
   cnv = prodc(abs(b - bold) .lt _bcriter) and (abs(devold - dev) lt _dcriter);
endo;
if iter gt _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
endif;
if bnames[1]; /* print results: */
   rreport(b,sqrt(diag(bcov)),bnames,yname,-1,-1,0,1);
   { dev0,score,inf } = mscore(xs,0,zeros(np,1),m,reps,ids);
   "Model degrees of freedom:";; format 5,0; np;
   "Deviance statistic and p for all coefficients:";;
   dev = dev0 - dev;
   format 11,3; dev;; format 11,5; cdfchic(dev,np);
   "Score statistic and p for all coefficients   :";;
   ss = score'invpd(inf)*score;
   format 11,3; ss;; format 11,5; cdfchic(ss,np);
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(b,bcov,dev,np);
endp;
