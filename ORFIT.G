/* This proc finds A-cells for 2x2 tables to give back
   the input odds ratios and margins using the quadratic formula.
    Input: r = odds ratios
           m = row-1 margins
           n = column-1 margins
           t = totals
    Output: (1,1)-cells that give back r given m,n,t */
/* Test program:
let r = 2 4; let m = 3 6; let n = 4 7; let t = 7 10;
r = r~(1/r);
a = orfit(r,m,n,t);
orf = a.*(t-m-n+a)./((m-a).*(n-a));
format /rds 6,2; r~orf;; a;; m~n~t; end; */
proc orfit(r,m,n,t);
local a2,bm,d,s;
/* a2 = 2*a, bm = -b, c = r.*m.*n: */
   a2 = 2*(r-1); bm = t + (r-1).*(m+n);
   d = sqrt(bm^2 - 2*a2.*r.*m.*n);
   s = d .le (a2.*min(m,n) - bm);
retp((bm + (2*s-1).*d)./(a2+(a2 .eq 0)) + (r .eq 1).*m.*n./t);
endp;
