/* This proc does Mantel-Haenszel analysis of person-time rates.
   Inputs: a = column vector of exposed cases
           b = col vector of unexposed cases
           t1 = col vector of exposed person-time
           t0 = col vector of unexposed person-time
           names = col vector of names of row, column, & stratifying variables
                   (set to zero if no output desired)
   Outputs: rdmh = Mantel-Haenszel rate difference -- default:
                   computed with the inefficient sparse-data weights
                   t1.*t0./(t1+t0); if _nullwt = 1, the null-efficient
                   weights t1.*t0./(a+b) will be used instead
                   (not valid if a+b averages less than 10)
            dse = Greenland-Robins SE for rdmh (Biometrics, 1985).
            rrmh = Mantel-Haenszel rate ratio.
            lse = Greenland-Robins SE for ln(rrmh) (Biometrics, 1985).
            chi = 1 df Mantel-Haenszel chi statistic (NOT CHI-SQUARED).
*/
proc (5) = mhrate(a,b,t1,t0,names);
local stab,m1,t,w,r,s,rdmh,dse,dd,ds,rs,ss,rrmh,vpt,lse,e,chi;
/* confidence coefficient: */ declare matrix _clevel = .95;
/* sparse-data indicator: */ declare matrix _nullwt = 0;
format /rds 8,4;
if ((a't0 eq 0) and (b't1 eq 0)) or (t1't0 eq 0);
   "MHRATE.G: No tables with nonzero diagonal or off-diagonal products.";
    retp(-1,-1,-1,-1,0);
endif;
/* eliminate zero-denominator tables: */
   stab = selif(a~b~t1~t0,(t1.*t0) .gt 0);
   a = stab[.,1]; b = stab[.,2]; t1 = stab[.,3]; t0 = stab[.,4];
/* margins: */ m1 = a+b; t = t1+t0;
if names[1]; print;
   ("MHRATE.G: Mantel-Haenszel analysis of the association of "
   $+names[1]$+" and "$+names[2]);;
   if rows(names) ge 3; format /rds 8,4;
      ","; " stratified by: ";; $names[3:rows(names)]';;
      "(";; format 1,0; rows(a);;"strata)";;
   endif; ".";
endif;
/* RD calculations: */
   if _nullwt; w = t1.*t0./m1; else; w = t1.*t0./t; endif;
   w = w/sumc(w); r = a./t1; s = b./t0; rdmh = w'(r-s);
/* RD SE: */ dse = sqrt((w.*w)'(r./t1 + s./t0));
/* eliminate strata with no cases: */ m1 = (a+b) .gt 0;
   a = selif(a,m1); b = selif(b,m1); t1 = selif(t1,m1); t0 = selif(t0,m1);
/* new margins: */ m1 = a+b; t = t1+t0;
/* null expectation & variance of sumc(a): */ e = m1.*t1./t; vpt = e.*t0./t;
/* rate ratio calculations: */ r = a.*t0./t; s = b.*t1./t;
   rs = sumc(r); ss = sumc(s);
   if ss; rrmh = rs/ss;
      if rs; /* ln(RR) SE calculations: */ lse = sqrt(sumc(vpt)./(rs.*ss));
         else; "WARNING: M-H rate ratio is zero."; lse = 0;
      endif;
      else; "WARNING: M-H rate ratio infinite or undefined.";
            rrmh = -1; lse = 0;
   endif;
/* M-H (score) statistic: */ chi = (sumc(a) - sumc(e))./sqrt(sumc(vpt));
if 0$+names[1]; /* print results: */ local z,rrcl,cnv,rrold,i,chih,df;
   z = cdfni(1-(1-_clevel)/2); z = 0~(-z)~z;
   "M-H estimate of uniform rate difference, and ";;
   format /rdn 1,0; 100*_clevel;;"% limits:";
   "          ";; format /rds 10,4; (rdmh + z*dse);;
      "; tanh limits:";; 
          dd = (1-rdmh)./(1+rdmh); ds = 2*z[2 3]*dse./(1-rdmh^2);
          (exp(ds)-dd)./(exp(ds)+dd);
   if rrmh gt 0;"M-H estimate of uniform rate ratio, and ";;
      format /rdn 1,0; 100*_clevel;; "% limits:";
      "          ";; format /rds 10,4; rrmh*exp(z*lse);
      /* find CMLE: */ rrcl = rrmh; rrold = 2*rrcl; cnv = 0; i = 0;
      do until i ge 20 or cnv; i = i + 1;
            cnv = abs(rrcl/rrold - 1) .lt .001;
            e = m1.*(1 - t0./(rrcl.*t1 + t0)); vpt = e.*(1 - e./m1);
            rrold = rrcl; rrcl = rrcl*exp(sumc(a-e)/sumc(vpt));
      endo;
      format /rdn 1,0; i;; "-step ";;
      if not cnv;"approximation to ";; endif;
      "CML estimate of uniform rate ratio, and ";;
      100*_clevel;; "% limits:";
      "          ";; format /rds 10,4; rrcl*exp(z/sqrt(sumc(vpt)));
   endif;
   "Mantel-Haenszel z-statistic & 2-sided p:";; format 8,3; chi;;
   format 8,5; cdfnt(chi);
   /* df for heterogeneity: */ df = rows(a)-1;
   format 8,0; "Numbers of exposed and unexposed cases:";; sumc(a~b)';
               "Minimum numbers of cases across strata:";; minc(a~b)';
   format 8,1; "Mean numbers of cases across strata   :";; meanc(a~b)';
   if rrmh gt 0 and df gt 0;
      "Approximate test for heterogeneity of the stratum-specific ratios";
      "(valid only for large strata) --";
      chih = (a - e)'((a - e)./vpt);
      " CML score statistic: ";; format 7,2; chih;; ", df = ";; format 3,0;
      df;;",  p =    ";; format 8,4; cdfchic(chih,df);
      if minc(b) .gt 0; "Range of stratum-specific ratios:";;
         format 9,2; minc(a.*t0./(b.*t1));; format 9,1; maxc(a.*t0./(b./t1));
      endif;
   endif;
endif;
format 9,3;
retp(rdmh,dse,rrmh,lse,chi);
endp;

