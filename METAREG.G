/* This program does log-linear weighted relative-risk regression,
   both fixed and random intercept ("random effects"),
   customized for meta-analyses of relative risks.
 Inputs: x = design matrix (meta-regressors),
         y = outcome vector (log relative risks),
         w = weight vector (0 if equal weights) or matrix
         const = 0 if no constant in regression,
         bnames = names of regressor (set to 0 if no printed output desired),
         yname = name of outcome (ignored if no printed output)
 Outputs: b = fixed-intercept coefficient estimates,
          bcov = inverse-information covariance matrix for b,
          y-x*b = residuals from fixed-intercept model,
          rss = weighted residual sum of squares for fixed effects,
          br = random-intercept coefficient estimates,
          bcovr = approximate covariance matrix for br,
          npa = estimated degrees of freedom for random-intercept model
                (rows(x)-npa is estimated residual df for model)
 */
proc (7) = metareg(x,y,w,const,bnames,yname);
   local t0,b,bcov,rss,rsq,df,v,br,bcovr,yp,sepy,t2,npa;
   t0 = date; if const; x = ones(rows(x),1)~x; endif;
   { b,bcov,rss,rsq,df } = wls(x,y,w,0,0,0,0);
   if rss gt df; /* fit random-intercept model: */
      if cols(w) eq 1; v = 1/w; else; v = invpd(w); endif;
      { br,bcovr,yp,sepy,t2,npa } = ebrr(y,v,x,0,0,0,0);
      else; br = b; bcovr = bcov; t2 = 0; npa = cols(x);
   endif;
   if 0$+bnames[1]; if const; bnames = (0$+"referent")|bnames; endif;
         "Proc METAREG.G: weighted log-linear regression";
         " with fixed and random intercept.";
         "No. of studies:";; format /rds 4,0; rows(x);
         "Fixed-effects results:";
         rreport(b,sqrt(diag(bcov)),bnames,yname,-1,-1,0,1);
         format 11,2; "HETEROGENEITY CHI-SQUARED, df, & p-value:";; rss;;
         format 9,0; df;; format 10,4; cdfchic(rss,df);
         if t2 gt 0; print;
            "Random-intercept results (assumes inverse-variance weights):";
             rreport(br,sqrt(diag(bcovr)),bnames,yname,-1,-1,-df,1);
             " Random-intercept SD estimate:";; format 10,3; sqrt(t2);
             else;
             "Random-intercept results not given because of low dispersion.";
          endif;
          "Total run time: ";; etstr(ethsec(t0,date));
   endif;
   retp(b,bcov,y-x*b,rss,br,bcovr,npa);
endp;
