/* Test program:  r = .2; n = 10;
let a = 4 1 3; eb = 0; call steinpr(a,n,r,eb,1);
let a = 4 1 3 5; eb = 1; call steinpr(a,n,r,eb,1); end; */
/* This proc computes the approximate Stein an EB estimators for risks & rates.
   Inputs: a = numerators
           n = denominators
           r = prior means (ignored if eb = 1)
           eb = 1 for EB, 0 for Stein
           c = 1 for count denominators, 0 for person-time
   Outputs: b = estimated posterior means
            w = estimated weights for a/n
*/
proc (2) = steinpr(a,n,r,eb,c);
declare matrix _maxit = 20;
local k,df,cd,cw,p,i,b,t2,v,d,u,w; print;
if sumc(n .eq 0); "INPUT ERROR IN STEINPR.G: ZERO DENOMINATOR"; end; endif;
k = rows(a); n = n.*ones(k,1);
if not eb; "PROC STEINPR.G: APPROXIMATE STEIN ESTIMATION";
   if k lt 3;
      "INPUT ERROR IN STEINPR.G: STEIN REQUIRES >2 OBSERVED NUMERATORS"; end;
   endif;
endif;
if eb; "PROC STEINPR.G: APPROXIMATE EMPIRICAL BAYES"; r = sumc(a)/sumc(n);
   if k lt 4;
      "INPUT ERROR IN STEINPR.G: EB REQUIRES >3 OBSERVED NUMERATORS"; end;
   endif;
endif;
r = r.*ones(k,1); df = k - eb; cd = k/df; cw = (df-2)/df; p = a./n;
i = 0; b = cw*p + (1-cw)*r; t2 = ((b - r)^2)./df;
do until i ge _maxit; i = i + 1;
   v = b.*(1-c.*b)./n;
   d = cd*(p-r)^2 - v;
   u = 1./(v + t2);
   t2 = (u'd)/sumc(u);
   w = 1-max(0,cw*u.*v);
   b = w.*p + (1-w).*r;
endo;
format /rds 6,0;
"Numerators       :";; a';
"Denominators     :";; n';
format /rds 6,2;
"Estimated weights:";; w';
"Observed ratios  :";; p';
"Prior means      :";; r';
"Posterior means  :";; b';
retp(b,w);
endp;
