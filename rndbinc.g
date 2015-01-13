/* test program: */
nsim = 1000; c = ones(nsim,1);
p=.1;
n=1000;
a = sumc(rndbinp(p,n,nsim) .eq k');
chi2 = ((o-e)^2)'(1/e);
chi2/(n+1);; cdfchic(chi2,n+1);; sumc(e~o)';; nsim;
e';o';
end;
proc rndbinp(p,n,nsim);
local k,pb,e,upb,lpb,u;
k = seqa(0,1,n+1);
pb = exp(lngamma(n+1)-lngamma(k+1)-lngamma(n-k+1)
          +k.*ln(p)+(n-k).*ln(1-p));
e=n.*pb;
upb = recserar(pb,pb[1,.],ones(1,cols(pb)));
lpb = 0|upb[1:n];
u=rndu(nsim,1);
retp( ((u .le upb').*(u .gt lpb'))*k );
endp;

/*
/* Generate a binomial variate 
 Code from C++ which is from R source that is from
 ACM Algorithm 678 BTPEC converted to C */

proc rndbinc(psave,n);
local p,q,np,r,g,qn,f,u,ix;

if n le 0 or psave le 0; retp(0); elseif psave ge 1; retp(n); endif;

if psave lt 0.5; p = psave; else; p = 1-psave; endif;
q = 1-p;
np = n*p;
r = p/q;
g = r*(n+1);

if np lt 30; 
  qn = q^n;
  do while 1;
    f = qn;
    u = rndu(1,1);
    ix = 0;
    do until ix gt 111; 
       if u lt f; goto finish; endif;
       ix = ix+1;
       u = u-f;
       f = f*(g/ix - r);
    endo;
  endo;
  goto finish; 
endif;

local m,c,ffm,fm,npq,p1,p2,p3,p4,xl,xll,xlr,xm,xr,x,al;
      ffm = np + p;
      m = trunc(ffm);
      fm = m;
      npq = np*q;
      p1 = trunc((2.195*sqrt(npq) - 4.6*q) + 0.5);
      xm = fm + 0.5;
      xl = xm - p1;
      xr = xm + p1;
      c = 0.134 + 20.5/(15.3 + fm);
      al = (ffm - xl)/(ffm - xl*p);
      xll = al*(1 + al/2);
      al = (xr - ffm)/(xr.*q);
      xlr = al*(1 + al/2);
      p2 = p1*(1 + c + c);
      p3 = p2 + c/xll;
      p4 = p3 + c/xlr;

  local i,v,k;
  do while 1;
    u = rndu(1,1)*p4; 
    v = rndu(1,1);
    if u le p1;
      ix = trunc(xm - p1*v + u);
      goto finish;
    endif;
    if u le p2;
      x = xl + (u - p1)/c;
      v = v*c + 1 - abs(xm - x)/p1;
      if v gt 1 or v le 0; continue; endif;
      ix = trunc(x);
    else;
      if u gt p3;
        ix = trunc(xr - ln(v)/xlr);
        if ix gt n; continue; endif;
        v = v*(u - p3) * xlr; 
      else; /* left tail */
        ix = trunc(xl + ln(v)/xll);
        if ix lt 0; continue; endif;
        v = v*(u - p2)*xll;
      endif;
    endif;

    k = abs(ix - m); 
    if k le 20 or k ge (npq/2 - 1);
      f = 1;
      if m lt ix; i = m+1; 
         do until i ge ix; i = i+1; f = f*(g/i - r); endo;
      elseif m ne ix; i = ix+1;
         do until i ge m; i = i+1; f = f*(g/i - r); endo;
      endif;
      if v le f; goto finish; endif;
      else;
      local amaxp, ynorm, alv;
      amaxp = (k/npq)*((k*(k/3 + 5/8)+ 1/6)/npq + 0.5);
      ynorm = -k*k/(2*npq);
      alv = ln(v);
      if alv lt ynorm - amaxp; goto finish; endif;
      if alv le ynorm + amaxp;
        local x1, f1, g1, z1, g2, x2, f2, z2;
        x1 = ix + 1;
        f1 = fm + 1;
        g1 = n - fm + 1;
        z1 = n - ix + 1;
        g2 = g1*g1;
        x2 = x1*x1;
        f2 = f1*f1;
        z2 = z1*z1;
        if alv le 
            xm*ln(f1/x1)+(n-m+0.5)*ln(g1/z1)+(ix-m)*ln(z1*p/(x1*q))
         +((13860 - (462 - (132 - (99 - 140/f2)/f2)/f2)/f2)/f1 
          +(13860 - (462 - (132 - (99 - 140/g2)/g2)/g2)/g2)/g1 
          +(13860 - (462 - (132 - (99 - 140/x2)/x2)/x2)/x2)/x1 
          +(13860 - (462 - (132 - (99 - 140/z2)/z2)/z2)/z2)/z1)/166320;
          goto finish; 
        endif;
      endif;
    endif;
  endo;

 finish:
  if psave gt 0.5; ix = n - ix; endif;
  retp(ix);
 endp;
*/