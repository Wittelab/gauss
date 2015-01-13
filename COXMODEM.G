/* This proc fits Cox models using the discrete-time Poisson approximation
   with a profile-likelihood algorithm based on the Nelson-Cox-Oakes
   approximation to the cumulative (discrete) hazard (p. 108 of Cox and Oakes).
   NOTES: 1. If n subjects are observed with k unique failure times, this
             proc will generate a k*n matrix, and several k*n row vectors.
          2. For time-dependent covariates, a separate record must be
             created for each covariate shift for each subject.
          3. Program does not calculate true deviances or loglikelihoods;
             those printed during iteration are for monitoring convergence
             and model comparisons only. They are not valid tests of fit.
          4. Baseline hazard and survival refer to distribution of failure
             times when all covariates are zero. Therefore, BE SURE YOUR
             COVARIATES ARE CODED SO THAT ZERO IS A MEANINGFUL VALUE FOR EACH!
   Inputs: x = matrix of covariates in model
           xs = matrix of excluded covariates to score test (0 if none)
           te = time of entry into observation (may be vector or a single time
                if everyone entered at the same time)
           tx = time of exit from observation (failure or censoring time)
           d = status at exit = 1 if failure, 0 if right-censored
           haz = 0 for no printout beyond coefficients & deviance,
                 1 for table of baseline survival,
                 2 for table & graphs of baseline survival
                 (set to 0 if case-control data!)
           rep = vector of repetion counts (0 if none)
           offset = vector of offsets (0 if no offsets)
           bnames = vector of fitted coefficient names
                  (set to zero if no printed output is desired)
           snames = vector of names for covariates to be score tested
                  (set to zero if no printed output is desired)
           yname = scalar name of outcome (may be 0 if no printed output).
 Outputs: b = coefficient estimates
          bcov = inverse-information covariance matrix
          dev = profile deviance for model
          sstat = score statistics for each column of xs
          ssall = score statistic for all columns of xs
Globals (may be set by user from calling program; otherwise defaults are used):
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _maxit = maximum number of iterations (default is 50)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .01)
      (WARNING: This algorithm can require a stringent deviance criterion)
*/
proc (5) = coxmodem(x,xs,te,tx,d,haz,rep,offset,bnames,snames,yname);
declare matrix _binit = 0; declare matrix _maxit = 50;
declare matrix _bcriter = .001; declare matrix _dcriter = .01;
   local t0,np,n,nones,nt,ft,start,atrisk,latrisk,tatrisk,dc,y;
   t0 = date;
   /* No. regressors: */ np = cols(x);
   if rank(x) lt np; "DESIGN MATRIX RANK DEFICIENT IN COXMODEM.G";
      retp(0,0,0,0,0); endif;
   /* no. of covariate levels: */ n = rows(x);
   /* constant vector of ones: */ nones = ones(n,1);
   if rows(rep) eq 1; y = d; rep = nones; else; y = rep.*d; endif;
   /* (y contains total number of deaths per record) */

      /* ******** BEGIN RISK-SET CONSTRUCTION ********* */
   /* If te is single number, set start and all entry times to te,
      where start is first entry time: */
      if rows(te) eq 1; start = te; te = te*nones;
         else; start = minc(te); endif; 
   /* failure times: */ ft = unique(selif(tx,d),1);
   /* no. of failure times: */ nt = rows(ft);
   /* Number of failures at each fail time: */ dc = sumc((tx .eq ft').*y);
   /* Matrix indicating presence of each subject at each failure time: */
      atrisk = (te .lt ft') .and (tx .ge ft');
   /* Vector giving index of last interval at risk for each subject: */
      latrisk = sumr(tx .ge ft') + (tx .lt ft[1]);
   /* Total at risk at each failure time: */ tatrisk = sumc(rep.*atrisk);
   bnames = 0$+bnames;
   if bnames[1]; print; /* output headers: */
      "PROC COXMODEM.G: Discrete Poisson approximation to Cox-model analysis.";
      if rows(bnames) ne cols(x);
        "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
      endif;
      format /rdn 6,0;
      "Number of subjects       : ";; sumc(rep);
      "Number of deaths         : ";; sumc(y);
      "Number of intervals      : ";; nt;
      "Number of regressors     : ";; np;
      format 10,3; "No data counted before time ";; start;;
                   " or after time ";; ft[nt]; print;
      "Preprocessing time (including assembly of risk sets): ";;
          etstr(ethsec(t0,date)); print;
    endif;

        /* ******** COMPUTE NULL STATISICS ********* */
   local hcrude,cumhaz,cscore,score,ssnull,llsat,llnull;
   /* Proportion failing at each ucutp: */ hcrude = dc./tatrisk;
   /* Initialize cumulative hazards for each observation at crude
      cumulative hazard (times repetition count): */
      cumhaz = rep.*(atrisk*hcrude);
   if bnames[1]; /* print null score statistic: */
      "Score test of all regressors   :";;
        cscore = ((x.*rep)'atrisk)'./tatrisk;
        score = (x'(y - cumhaz));
        ssnull = score'invpd((x.*cumhaz)'x - (cscore.*dc)'cscore)*score;
      format 11,3; ssnull;;
      format 7,5; ", p = ";; cdfchic(ssnull,np);;
      format /rdn 4,0; ", df = ";; np; print;
   endif;
   /* Saturated (= crude) Poisson loglikelihood: */
      llsat = dc'ln(dc) - sumc(dc) - sumc(lnfact(dc));
   /* Null Poisson loglikelihood: */
      llnull = y'ln(hcrude[latrisk]) - sumc(cumhaz);

        /* ******* START MODEL FITTING ******* */
   local dev,devold,b,bold,cnv,iter,trhx,trhf,hdisc,eta,mu,inf,bcov;
   /* Temporarily augment x with constant: */ x = nones~x;
   /* Temporarily augment no. of fitted parameters: */ np = cols(x);
   /* Initialize loop for model fitting: */
      if (rows(_binit) eq 1) and (_binit eq 0); b = zeros(np,1);
         else; b = _binit; endif;
      dev = 0; cnv = 0; iter = 0;
   do until (iter ge _maxit) or cnv; iter = iter + 1;
        /* save old values: */ bold = b; devold = dev;
        /* linear predictor: */ eta = offset + x*b;
        /* estimate the baseline hazard given beta: */
           /* total relative hazard at each covariate (x) level: */
              trhx = rep.*exp(eta);
           /* total relative hazard (adj. no. at risk) at each failure time: */
              trhf = (trhx'atrisk)';
           /* estimated discrete baseline hazards: */ hdisc = dc./trhf;
           /* estimated cumulative baseline hazards: */ cumhaz = atrisk*hdisc;
        /* expected case count for each covariate level: */ mu = cumhaz.*trhx;
        dev = -2*(y'(eta + ln(hdisc[latrisk])) - sumc(mu) - llsat);
        if bnames[1]; "Deviance at iteration ";;
           format /rdn 4,0; iter;; ":";; format 11,3; dev; endif;
        /* estimate the betas given the baseline hazards: */
           /* information matrix: */ inf = moment(sqrt(mu).*x,0);
           /* inverse information: */ trap 1; bcov = invpd(inf);
           if scalerr(bcov); "SINGULARITY IN COXMODEM.G";
              retp(0,0,0,0,0); endif; trap 0;
           /* Score vector: */ score = x'(y - mu);
           /* Newton step: */ b = b + bcov*score;
        /* test for convergence: */
        cnv = prodc(abs(b-bold) .< _bcriter) and (abs(devold-dev) < _dcriter);
   endo;
   if iter ge _maxit;
      "WARNING: No convergence after ";; iter;; " iterations"; endif;

             /* ******* POSTFIT STATISTICS ******* */
   local devnull,cscs,sstat,ssall,s0;
   /* rescale estimated discrete hazard: */ hdisc = exp(b[1])*hdisc;
   /* remove intercept: */
      x = x[.,2:np]; inf = inf[2:np,2:np]; score = score[2:np];
      b = b[2:np]; np = cols(x);
   if bnames[1];
      print; "Deviance test of all regressors:";;
      devnull =  2*(llsat - llnull) - dev;
      format 11,3; devnull;;
      format 7,5; ", p = ";; cdfchic(devnull,np);;
      format /rdn 4,0; ", df = ";; np;
   endif;
   /* Compute covariance matrix for b: */
      /* Total conditional score for each interval: */
         cscore = ((x.*trhx)'atrisk)'./trhf;
      /* Weighted sscp matrix of conditional scores: */
         cscs = (cscore.*dc)'cscore;
      /* inf is the total information;
         cscs is the information used in estimation of the hazard function;
         hence the information left for b is inf - cscs: */
     /* Covariance matrix of b: */ bcov = invpd(inf - cscs);
   sstat = -1; ssall = -1;
   if rows(xs) .eq rows(x); /* compute score tests on columns of xs: */
      if snames[1]; print; " Score statistics for adding"; endif;
      /* Preweight cscore by dc: */ cscore = cscore.*dc;
      sstat = {}; iter = 0;
      do until iter ge cols(xs); iter = iter + 1;
         sstat =
   sstat|coxscore(x,xs[.,iter],y,dc,mu,trhx,atrisk,trhf,score,inf,cscore,cscs);
         if snames[1];
            "  ";; format /ldn 10,0; $snames[iter];;
            ":";; format /rdn 9,3; sstat[iter];;", p = ";;
            format 7,5; cdfchic(sstat[iter],1);;
            format 5,0; ", df = ";; cols(xs);
         endif;
      endo;
      if cols(xs) gt 1; /* do summary score test: */
         "  ";; format /ldn 10,0; "all above :";;
         ssall = coxscore(x,xs,y,dc,mu,trhx,atrisk,trhf,score,inf,cscore,cscs);
         format /rdn 9,3; ssall;;", p = ";; format 7,5;
         cdfchic(ssall,cols(xs));; ", df = ";; format 5,0; cols(xs);
      endif;
   endif;

             /* ******* FINAL REPORT ******* */
if bnames[1]; /* report results: */ print;
  rreport(b,sqrt(diag(bcov)),bnames,yname,-1,-1,0,1); print;
  "Total run time to this point: ";; etstr(ethsec(t0,date)); print;
  if haz; /* report baseline hazards and survival dist: */
    s0 = exp(-cumsumc(hdisc));
    "Failure times and number of failures, number at risk,";
    " proportion failing, estimated baseline discrete hazard,";
    " and estimated baseline survival probability at that time:";
    format 9,3; ft~dc~tatrisk~hcrude~hdisc~s0;
    if haz gt 1; /* plot fitted hazards and survival dist: */
      library pgraph;
      "Press enter to see plot of estimated baseline survival distribution.";
      wait; _pltype = 4;
 /* This is a trick to get Gauss to graph a step function: */
 xy(start|(ft.*.ones(2,1)),1|reshape((1|s0[1:(nt-1)])~s0,2*nt,1));
      "Press enter to see connected plot of estimated baseline survival dist.";
      wait; xy(start|ft,1|s0);
    endif;
  endif;
endif;
retp(b,bcov,dev,sstat,ssall);
endp;

proc coxscore(x,xs,y,dc,mu,trhx,atrisk,trhf,score,inf,cscore,cscs);
local xsmu,infs,cxscore,cscxs,sstat;
         xsmu = xs.*mu;
         /* Weighted sscp of x with xs: */ infs = x'xsmu;
         /* Augment total wtd sscp matrix with infs and xswxs: */
         inf = (inf~infs)|(infs'~(xsmu'xs));
         /* Augment x-scores with xs-score: */ score = score|(xs'(y-mu));
         /* Conditional xs-scores: */
            cxscore = ((xs.*trhx)'atrisk)'./trhf;
         /* Weighted sscp of cond x-scores with cond xs-scores: */
            cscxs = cscore'cxscore;
         /* Augment cscs with cscxs and cxscore'cxscore: */
            cscs = (cscs~cscxs)|(cscxs'~((cxscore.*dc)'cxscore));
         /* Score stat for adding xs: */
            sstat = score'invpd(inf - cscs)*score;
retp(sstat);
endp;
