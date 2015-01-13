/* This proc does the Robins test on br and bp.
   Inputs:
    rind = treatment-arm indicator
    tc = withdrawal time (known censoring time)
    tx = min(tc,td) (death or death-censoring time)
    tcpd = min(tc,tp,td) = min(tp,tx) (prophylaxis or proph-censoring time)
    br = test value for brt
    bp = test value for bpt
    prin = 1 if printed output desired, 0 if not
   Intermediates:
    tz = transformed observable time (z-dagger)
       = t00|br,bp if less than all blip-downs of tc
    blip*tc = censoring time for analysis of tz
    tu = minimum blip-down time
    uind = tz .lt blip*tc (indicates tz is not censored by blip*tc,
           in which case tz = t00|br,bp)
   Outputs:
    lrs = total logrank score
    lrv = total logrank variance
    pws = total Prentice-Wilcoxon score
    pwv = total PW variance
    bwt = bivariate test statistic based on logrank and PW scores
*/

    proc (5) = aftg(rind,tc,tx,tcpd,br,bp,prin);
    local blip,reff,peff,tz,tu,uind,lrs,lrv,lrt,pws,pwv,bwt;
    /* Logrank test procedures: */
    blip = exp(minr(0~bp~br~(br+bp)));
    reff = exp(br*rind);
    peff = exp(bp);
    tz = reff.*(tcpd + peff*(tx - tcpd));
    tu = minr((blip*tc)~tz);
    uind = tz .lt blip*tc;
    if prin;
       format 9,3; "For beta-r, beta-p = ";; br~bp;;"   --";
    endif;
    { lrs,pws,lrv,pwv,bwt } = logrank(rind,tu,uind,prin);
retp(lrs,pws,lrv,pwv,bwt);
endp;