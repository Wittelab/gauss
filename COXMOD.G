/* This proc fits Cox models using exponential regression; it optionally
   tabulates and graphs the corresponding baseline survival curve.
   RUNS VERY SLOWLY COMPARED TO COXMODEM.G
   NOTES: 1. If n subjects are observed an average of k intervals each, this
             proc will generate k*n row working data matrices.
          2. For time-dependent covariates, a separate subject record must be
             created for each covariate shift for each subject, with entry
             time at the covariate shift.
          3. Baseline hazard and survival refer to distribution of failure
             times when all covariates are zero. Therefore, BE SURE YOUR
             COVARIATES ARE CODED SO THAT ZERO IS A MEANINGFUL VALUE FOR EACH!
   Inputs: x = covariate matrix
           te = time of entry under observation
           tx = time of exit from observation
           d = status at exit = 1 if failure, 0 if right-censored
           lcutp = start of first interval or starts of all intervals
           ucutp = ends of all intervals or end of last interval (set to 0
                   if you want an interval between each event time tx, with
                   first interval starting at earliest entry time)
           tpow = highest power of time to model in hazard -
                  if ge 0, polynomial regression will be used;
                  if lt 0, time will be treated categorically, which
                           produces nonparametric estimation of baseline hazard.
           haz = 1 if estimation of baseline hazard and survival curve
                 desired, 0 otherwise.
           a = vector of repetion counts (0 if none),
           offset = vector of offsets (0 if no offsets),
           bnames = vector of coefficient names
                  (set to zero if no printed output is desired),
           yname = scalar name of outcome (may be 0 if no printed output).
 Outputs: b = coefficient estimates,
          bcov = inverse-information covariance estimate,
          bcovr = specification-robust covariance estimate,
          dev = pseudo-residual deviance for model,
          df = pseudo-residual degrees of freedom (for model comparison only).
Globals (may be set by user from calling program; otherwise defaults are used)
      -- these are used in expreg.g, which is called in this procedure:
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _maxit = maximum number of iterations (default is 30)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
*/
proc (5) = coxmod(x,te,tx,d,lcutp,ucutp,tpow,haz,a,offset,bnames,yname);
local t0,n,nones,events,nin,lnth,nin1;
t0 = date;
   n = rows(x); nones = ones(n,1);
   /* Indicates if intervals are to be constructed from event times: */
      events = sumc(ucutp) eq 0;
if events;
   if tpow ge 0; /* Use all event times for polynomial regression: */
         ucutp = unique(tx,1);
      else; /* Use only failure times for categorical time: */
         ucutp = unique(selif(tx,d),1);
   endif;
   /* Lower cutpoints: */ lcutp = minc(te)|ucutp[1:(rows(ucutp)-1)];
endif;
/* Interval lengths: */ lnth = ucutp - lcutp;
/* No. of intervals: */ nin = rows(ucutp); nin1 = ones(nin,1);

/* Construct records for each person-interval: */
   local nucutp,txn,ppt,dn,an,m1,offn;
   /* Upper cutpoints for each person-interval: */ nucutp = nones.*.ucutp;
   /* End of observation for each person-interval: */ txn = tx.*.nin1;
   /* If everyone starts at same time, make a vector out of this time: */
      if rows(te) eq 1; te = te*nones; endif;
   /* Calculate observed time for each person-interval: */
   if events and (tpow lt 0); /* Count time only in intervals with failures: */
      ppt = (((txn .ge nucutp) .and (te.*.nin1 .lt nucutp)).*(nones.*.lnth));
   else; /* Count time for each person-interval: */
      ppt = pospart(minr(nucutp~txn) - maxr((nones.*.lcutp)~(te.*.nin1)));
   endif;
   /* Death indicators for each person-interval: */
      dn = (txn .eq nucutp) .and (d.*.nin1);
   /* Repetition counts (an) for each person interval
      and total no. deaths (m1) in each interval: */
   if sumc(a); an = selif(a.*.nin1,ppt);  m1 = sumc(a.*reshape(dn,n,nin));
               else; an = 0; m1 = sumc(reshape(dn,n,nin)); endif;
   if sumc(offset); offn = selif(offset.*.nin1,ppt); else; offn = 0; endif;

/* Fit model: */
   local tvar,tnames,b,bcov,bcovr,dev,rss,df;
   if sumc(a) eq 0; a = nones; endif;
   if tpow ge 0; /* Do polynomial Poisson regression: */
       /* First, construct time variable(s); time is centered: */
       tvar = (ucutp - meanc(ucutp))^seqa(0,1,tpow+1)';
     else;
       /* Categorical time (WARNING: May overflow memory if many cutpoints): */
       tvar = eye(nin);
   endif;
   bnames = 0$+bnames;
   if bnames[1];
      "PROC COXMOD: Cox-Poisson failure-time analysis";
      if rows(bnames) ne cols(x);
        "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
      endif;
endif;
   { b,bcov,bcovr,dev,rss,df } =
         expreg(selif(nones.*.tvar,ppt)~selif(x.*.nin1,ppt),
                selif(dn,ppt),selif(ppt,ppt),
                an,0,offn,0,0);

if bnames[1]; /* Report results: */ local h0,se;
   /* Baseline hazard: */ h0 = exp(tvar*b[1:cols(tvar)]);
   se = sqrt(diag(bcov));
   if tpow ge 0; /* Include time coefficients in report: */
         " with polynomial rate model of order ";; format /rdn 2,0; tpow;;".";
         let tnames = "Const" "t" "t^2" "t^3" "t^4" "t^5" "t^6";
         bnames = tnames[1:(tpow+1)]|bnames;
      else; /* Exclude time coefficients from report: */
         " with categorical time (nonparametric step-function hazard).";
         b = b[(nin+1):rows(b)]; se = se[(nin+1):rows(se)];
   endif;
   format /rdn 6,0;
   "Number of subjects       : ";; sumc(a);
   "Number of deaths         : ";; a'd;
   "Number of intervals      : ";; nin;
   format 10,3; "No data counted before time ";; lcutp[1];;
                " or after time ";; ucutp[nin];
   rreport(b,se,bnames,yname,dev,-1,df,1); print;
   "Total run time: ";; etstr(ethsec(t0,date));

 if haz; /* Compute, tabulate, and plot survival estimates: */
   library pgraph; local s0,pat,nat; format 9,3;
   /* Baseline survival prob, and no. at risk in each interval: */
      s0 = exp(-cumsumc(h0.*lnth));
   /* total person-time at risk in each interval: */
      ppt = reshape(ppt,n,nin);
      pat = sumc(a.*ppt);
      nat = sumc(a.*(ppt .gt 0));
   "Interval boundaries, failures, no. at risk at start, person-time at risk,";
   " crude failure rates (failures/person-time at risk),";
   " estimated baseline failure rates, and estimated survival probabilities:";
   " (intervals with no failures are omitted if event boundaries and";
   "  categorical time are used simultaneously, for then the crude and";
   "  estimated baseline hazards will be zero in intervals with no failures)";
   lcutp~ucutp~m1~nat~pat~(m1./pat)~h0~s0;
   "Press enter to see plot of estimated baseline survival distribution.";
   wait; _pltype = 4;
   xy(lcutp[1]|(ucutp.*.ones(2,1)),1|reshape((1|s0[1:(nin-1)])~s0,2*nin,1));
   "Press enter to see connected plot of estimated baseline survival dist.";
   wait;
   xy(lcutp[1]|ucutp,1|s0);
 endif;

print; endif;
retp(b,bcov,bcovr,dev,df);
endp;



