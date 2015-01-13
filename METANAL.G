/* This proc does basic unstratfied meta-analysis -
   see METAREG.G for meta-regression.
   Inputs: b = study-specific estimates (e.g., ln(rr))
           w = study-specific weights (usually inverse variances)
           f = 1 for funnel plot, 0 not
           xname = study covariate (0 if no printed output desired)
           yname = study outcome (0 if no printed output desired)
   Outputs: bf,1/ws = fixed-effects summary & its variance estimate
            brs,vrs,brh,vrh,brm,vrm = random-effects summaries & their
            variance estimates from DerSimonian-Laird, Hardy-Thompson, & Morris
            t2s,t2h,t2m = estimates of variance of RE distribution
   Global: _nopub = 1 if no regression tests of publication bias desired
*/
proc (11) = metanal(b,w,f,xname,yname);
declare matrix _nopub = 0;
local n,df,ws,bf,xh,t2s,wr,vrs,brs,i,c,t2old,brh,wrh,vrh,t2h,
      brm,vrm,yp,sepy,t2m,npa,rad,t,rt,rw,mt,mw,x,z;
n = rows(b); df = n-1; ws = sumc(w);
/* fixed-effects estimate: */ bf = (w'b)/ws;
/* het chi-squared: */ xh = w'((b - bf)^2);
/* DerSimonian-Laird RE variance estimate: */
   t2s = max(ws*(xh - df)/(ws*ws - w'w),0);
   if t2s le 0; "DerSimonian-Laird yields underdispersion."; endif;
/* RE weights: */ wr = 1/(1/w + t2s); vrs = 1/sumc(wr);
/* DerSimonian-Laird RE estimate: */ brs = (wr'b)*vrs;
/* Hardy-Thompson iterative RE estimates (ML for normal b): */
   i = 0; c = n/df; t2old = -1; t2h = t2s; wrh = wr; brh = brs;
   do until i ge 50 or abs(t2h-t2old) lt .0001; i = i + 1;
      wrh = wrh.*wrh; t2old = t2h;
      t2h = max((wrh'(c*(b-brh)^2 - 1/w))/sumc(wrh),0);
      wrh = 1/(1/w + t2h); brh = (wrh'b)/sumc(wrh);
   endo;
   vrh = 1/sumc(wrh);
   if t2h le 0; "Hardy-Thompson yields underdispersion."; endif;
/* Carl Morris 1983 JASA estimates: */
   { brm,vrm,yp,sepy,t2m,npa } = ebrr(b,1/w,ones(n,1),0,0,0,0);
   if t2m le 0; "Carl Morris yields underdispersion."; t2m = 0; endif;
if 0$+xname; ("Unstratified meta-analysis of "$+ftocv(n,1,0)$+" studies of "
   $+xname$+" and "$+yname);
   "Chi-squared for heterogeneity:";; format /rds 9,2; xh;
   ""$+ftocv(n-1,1,0)$+" degrees of freedom, P =";;
   format 8,5; cdfchic(xh,df); format 9,2;
   let rad[1,3] = 0 -1 1;
   "Fixed effects summary & 95% CL :";; exp(bf + 1.96*rad/sqrt(ws));;
    format 8,5; ", P =";; cdfnt(bf*sqrt(ws)); format 9,2;
   print;
   if t2s gt 0;
      "Random effects summary & 95% CL based on Dersimonian-Laird & t dist:";
      "  ";; exp(brs + invt(.95,df)*rad*sqrt(vrs));;
      format 8,5; ", P =";; 2*cdftc(abs(brs/sqrt(vrs)),df); format 9,2;
      "  95% interval for random-effects distribution:";;
      exp(brs + 1.96*rad[2 3]*sqrt(t2s));
      print;
   endif;
   if t2h gt 0;
      "Random effects summary & 95% CL based on Hardy-Thompson & t dist:";
      "  ";; exp(brh + invt(.95,df)*rad*sqrt(vrh));;
      format 8,5; ", P =";; 2*cdftc(abs(brs/sqrt(vrh)),df); format 9,2;
      "  95% interval for random-effects distribution:";;
      exp(brh + 1.96*rad[2 3]*sqrt(t2h));
      print;
   endif;
   if t2s gt 0;
      "Random effects summary & 95% CL based on Morris & t dist:";
      "  ";; exp(brm + invt(.95,df)*rad*sqrt(vrm));;
      format 8,5; ", P =";; 2*cdftc(abs(brm/sqrt(vrm)),df); format 9,2;
      "  95% interval for random-effects distribution:";;
      exp(brm + 1.96*rad[2 3]*sqrt(t2m));
      print;
   endif;
   "Begg-Mazumdar test for publication bias: P = ";;
    t = (b-bf)./sqrt(1/w - 1/ws);
    /* rankings: */ rt = rankindx(t,1); rw = rankindx(w,1);
    /* rank comparisons: */ mt = rt .gt rt'; mw = rw .gt rw';
    /* no. concordant rankings: */ x = sumall(upmat1(mt .eq mw))-n;
    z = (2*x - n*df/2)/sqrt(n*df*(2*n+5)/18);
    format 8,5; cdfnt(z); format 9,2;
   if not _nopub;
      print; "Regression tests of asymmetry of funnel plot";
      " (slope should be zero if symmetric):";
      call wls(w,b,1,0,1,"weight",yname);
      call wls(w,t,1,0,1,"weight","t-score");
   endif;
endif;
if f; library pgraph; _pdate =""; fonts("microb");
   let _pstype = 9; let _pcolor = 3; let _psymsiz = 3.5; let _plctrl = -1;
   let _pline[1,9] = 1 4 0 0 0 200 0 7 9; _pline[1,6] = maxc(w);
   xlabel("Funnel plot of "$+ftocv(n,1,0)$+" estimates");
   ylabel("study weight"); xy(b,w);
endif;
retp(bf,1/ws,brs,vrs,brh,vrh,brm,vrm,t2s,t2h,t2m);
endp;
