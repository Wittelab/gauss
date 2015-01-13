/* This proc does unweighted (1 df) and Prentice-Wilcoxon (1 df),
   and returns the test scores and variances.
   Inputs: treat = treatment or exposure variable
           tx = min(failure time, censoring time)
           dind = 1 if tx is failure time, 0 if censoring time
           prin = 1 if printed output desired, 0 if not
   Intermediates:
    txt = tx among the treated
    utd = vector of unique failure times
    nat = no. at risk at utd (no. with tx ge utd)
    mts = mean of treatment scores at utd
    mat = no. failures at utd
    s = vector of logrank scores
    v = vector of logrank variances
    kmu = Kaplan-Meier survival estimates at utd for pooled data
    lrt = logrank Z statistic
    pwt = PW Z statistic
   Outputs:
    lrs = total logrank score
    lrv = total logrank variance
    pws = total Prentice-Wilcoxon score
    pwv = total PW variance
*/
    proc (4) = logrankp(treat,tx,dind,prin);
    local n,ind,utd,nat,mat,mts,s,v,lrs,lrv,lrt,pws,pwv,pwt,cov,lrcov,kmu;
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
    lrt = lrs/sqrt(lrv);
    kmu = cumprodc(1 - (mat./nat));
    pws = kmu's;
    cov = kmu.*v; lrcov = sumc(cov);
    pwv = kmu'cov;
    pwt = pws/sqrt(pwv);
    if prin;
      "No. 'failures' & no. unique failure times :";; 
       format 7,0; sumc(dind);; rows(utd);
      "Estimated logrank and Wilcoxon correlation:";; 
       format 9,3; lrcov/sqrt(lrv*pwv);
      "Logrank  Z statistic and P:";; format 9,3; lrt;; format 8,5; cdfnt(lrt);
      "Wilcoxon Z statistic and P:";; format 9,3; pwt;; format 8,5; cdfnt(pwt);
       format 9,3; print;
     endif;
retp(lrs,lrv,pws,pwv);
endp;