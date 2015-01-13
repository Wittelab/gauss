/* This proc does unweighted (1 df), Prentice-Wilcoxon (1 df),
   and Robins (2 df) logrank trend tests on input data,
   and returns the test scores and variances.
   Inputs: treat = treatment or exposure variable
           tx = min(failure time, censoring time)
           dind = 1 if tx is failure time, 0 if censoring time
           prin = 1 if printed output desired, 0 if not
   Intermediates:
    txt = tx among the treated
    utd = unique(txz,1) = vector of unique failure times
    nat = no. at risk at utd (no. with tx ge utd)
    mts = mean of treatment scores at utd
    mat = no. failures at utd
    s = vector of logrank scores
    v = vector of logrank variances
    kmu = Kaplan-Meier survival estimates at utd for pooled data
    cov = covariance of lrs and pws
    lrt = logrank test statistic
    pwt = PW test statistic
   Outputs:
    lrs = total logrank score
    lrv = total logrank variance
    pws = total Prentice-Wilcoxon score
    pwv = total PW variance
    bwt = bivariate test statistic based on logrank and PW scores
*/
    proc (5) = logrank(treat,tx,dind,prin);
    local n,ind,utd,nat,mat,mts,s,v,
          lrs,lrv,lrt,pws,pwv,pwt,bws,bwv,bwt,kmu,cov,lrcov,df;
    n = rows(tx);
/* Sort treat, tx, and dind on tx: */
    ind = sortind(tx); treat = treat[ind]; tx = tx[ind]; dind = dind[ind];
    utd = unique(selif(tx,dind),1);
    nat = countge(tx,utd,0);
    mts = countge(tx,utd,treat)./nat;
    mat = countwts(tx,utd,dind);
    s = countwts(tx,utd,dind.*treat) - mat.*mts;
    v = countge(tx,utd,treat.*treat)./nat - mts.*mts;
    v = v.*(nat - mat).*mat./(nat - (nat .gt 1));
    lrs = sumc(s); lrv = sumc(v);
    lrt = lrs*lrs/lrv;
    kmu = cumprodc(1 - (mat./nat));
    pws = kmu's;
    cov = kmu.*v; lrcov = sumc(cov);
    pwv = kmu'cov;
    pwt = pws*pws/pwv;
    bws = lrs|pws;
    bwv = invpd((lrv|lrcov)~(lrcov|pwv));
    bwt = bws'bwv*bws;
    if prin; print; format 9,3;
      "No. 'failures' & no. unique 'failure times':";; sumc(dind);; rows(utd);
       format /rdn 9,3;" Logrank and P-W scores:";; lrs~pws;
       " Logrank and P-W covariance matrix and inverse:";
         ((lrv|lrcov)~(lrcov|pwv)); bwv;
  " Logrank, Prentice-Wilcoxon, and Robins test statistics, df, and p values:";
         format 9,3; lrt~pwt~bwt; let df[1,3] = 1 1 2; format 9,0; df;
         format 9,5; cdfchic(lrt~pwt~bwt,df); format 9,3; print;
     endif;
retp(lrs,pws,lrv,pwv,bwt);
endp;

