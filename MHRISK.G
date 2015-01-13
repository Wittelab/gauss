/* This proc does Mantel-Haenszel analysis of average risks
   (incidence proportions) in cohort data.
   Inputs: a = column vector of exposed cases
           b = col vector of unexposed cases
           n1 = col vector of exposed totals
           n0 = col vector of unexposed totals
           names = col vector of names of row, column, & stratifying variables
                   (set to zero if no output desired)
   Outputs: rdmh = Mantel-Haenszel risk difference -- default:
                   computed with the inefficient sparse-data weights
                   n1.*n0./(n1+n0); if _nullwt = 1, the more efficient
                   weights n1.*n0./(a+b) will be used instead
                   (not valid if a+b averages less than 10)
            dse = binomial SE for rdmh -
                  computed only if n1 and n0 are always 2 or more
            rrmh = Mantel-Haenszel risk ratio -- default:
                   computed with the inefficient sparse-data weights
                   n1.*n0./(n1+n0); if _nullwt = 1, the more efficient
                   Tarone weights n1.*n0./(c+d) will be used instead
                   (not valid if c+d averages less than 10)
            lse = Greenland-Robins SE for sparse-data ln(rrmh)
                  (Biometrics, 1985).
            rree, lsee = estimating-equation risk ratio estimate and estimated SE, based on iterative refinement of MH weights
            chi = 1 df Mantel-Haenszel chi statistic (NOT CHI-SQUARED).
*/
proc (7) = mhrisk(a,b,n1,n0,names);
local stab,m1,t,w,r,s,i,rdmh,dse,ds,dd,rs,ss,rrmh,vpt,lse,rree,rrold,lsee,e,chi;
/* confidence coefficient: */ declare matrix _clevel = .95;
/* sparse-data indicator: */ declare matrix _nullwt = 0;
format /rds 8,4;
if ((a'n0 eq 0) and (b'n1 eq 0)) or (n1'n0 eq 0);
   "MHRISK.G: No tables with nonzero diagonal or off-diagonal products.";
   retp(-1,-1,-1,-1,-1,-1,0);
endif;
/* eliminate zero-denominator tables: */
   if sumc(n1.*n0) eq 0; "MHRISK.G: Zero denominators in all tables!";
      retp(0,-1,-1,-1,0); endif;
   stab = selif(a~b~n1~n0,(n1.*n0) .gt 0);
   a = stab[.,1]; b = stab[.,2]; n1 = stab[.,3]; n0 = stab[.,4];
/* margins: */ m1 = a+b; t = n1+n0;
/* RD calculations: */
   if _nullwt; w = n1.*n0./m1; else; w = n1.*n0./t; endif;
   w = w/sumc(w); r = a./n1; s = b./n0; rdmh = w'(r-s);
/* RD SE: */
   if sumc((n1|n0) .lt 2); dse = -1; else;
      dse = sqrt((w.*w)'(r.*(1-r)./(n1-1) + s.*(1-s)./(n0-1))); endif;

/* eliminate strata with no cases: */ m1 = (a+b) .gt 0;
   a = selif(a,m1); b = selif(b,m1); n1 = selif(n1,m1); n0 = selif(n0,m1);
/* new margins: */ m1 = a+b; t = n1+n0;
/* null expectations of a: */ e = m1.*n1./t; vpt = e.*n0./t;
/* M-H (score) statistic: */  chi = sumc(a-e)./sqrt(vpt'((t-m1)./(t-1)));
/* RR calculations: */
   if _nullwt; r = a.*n0./(t-m1); s = b.*n1./(t-m1);
      else; r = a.*n0./t; s = b.*n1./t; endif;
   rs = sumc(r); ss = sumc(s);
   if ss; rrmh = rs/ss; rree = rrmh;
      if rs; /* ln(RR) SE calculations: */
         lse = sqrt(sumc(vpt - a.*b./t)/(rs*ss));
         i = 0;
         do until i ge 30; i = i + 1;  rrold = rree;
              w = 1/(rree*n1+n0);
              rs =  w'(a.*n0); ss = w'(b.*n1); rree = rs/ss;
              if abs(ln(rree/rrold)) lt 0.0001; continue; endif;
         endo;
         lsee = sqrt(((w.*w)'(m1.*n1.*n0 - a.*b.*t))/(rs*ss)); 
         else; "WARNING: M-H and EE risk ratios are zero"; lse = 0; lsee = 0;
      endif;
      else; "WARNING: M-H and EE risk ratios are  infinite or undefined";
            rrmh = -1; rree = -1; lse = 0; lsee = 0;
   endif;
if 0$+names[1]; /* print results: */ local z,chih,df;
   z = invnorm(1-(1-_clevel)/2); z = 0~(-z)~z;
   ("MHRISK.G: Mantel-Haenszel analysis of the association of "
   $+names[1]$+" and "$+names[2]);;
   if rows(names) ge 3; format /rds 8,4;
      ","; " stratified by: ";; $names[3:rows(names)]';;
      "(";; format 1,0; rows(a);;"strata)";;
   endif; ".";
   if _nullwt; " -- Using null weights."; endif;
   "Estimated uniform risk difference";;
   if dse le 0; ":";; format 10,5; rdmh;
      else; ", and ";; format /rdn 1,0; 100*_clevel;;"% limits:";
      "          ";; format /rds 10,5; (rdmh + z*dse);;
      "; tanh limits:";; 
          dd = (1-rdmh)./(1+rdmh); ds = 2*z[2 3]*dse./(1-rdmh^2);
          (exp(ds)-dd)./(exp(ds)+dd);
   endif;
   "MH and EE estimates of uniform risk ratio, and ";; format /rdn 1,0;
   100*_clevel;; "% limits"; format /rds 10,4;
   "          ";; rrmh*exp(z*lse); "          ";; rree*exp(z*lsee); format 8,2;
   "Mantel-Haenszel z-statistic & 2-sided p:";; chi;; format 8,5; cdfnt(chi);
   if _nullwt and prodc((a|b) .ne 0);
      "Approximate test for heterogeneity of the stratum-specific ratios";
      "(valid only for large strata and summary ratio not too far from 1):";
      chih = (ln(a.*n0./(b.*n1)) - ln(rree)); chih = chih'(chih./(1/a + 1/b));
      " M-H Wald statistic: ";; format 7,2; chih;; ", df = ";; format 1,0;
      df = rows(a)-1; df;;",  p =    ";; format 8,5; cdfchic(chih,df);
   endif;
endif;
format 9,3;
retp(rdmh,dse,rrmh,lse,rree,lsee,chi);
endp;
