/* Test program: */
"From Breslow & Day 1980, p. 126:";
format /rdn 10,4;
p = phyperg(3,5,7,10|1); p;
"lower p values for OR = 10, 1: ";; sumr(p[.,2 3])';
print; "From Logexact ad of Summer 1993:";
p = phyperg(36,5,1086,24.5|29.1); p;
print; "Recap:"; p = phyperg(3|36,5|5,7|1086,25); p;
end;
/* This proc computes the noncentral hypergeometric distribution for the
   A cells in 2-by-2 tables with totals m and n in the A row and column,
   grand totals t, and odds ratios r.
   The returned distributions are the rows of the matrix p; each row gives
   probabilities for A = 0 to A = gmax, where gmax is the largest possible
   value of A over all the tables.
   Element k,u of p is thus the probability that A = u-1 for table k.
*/
proc phyperg(m,n,t,r);
local amin,amax,gmax,k,cc,c,e,s;
amin = max(0,m-t+n);
amax = min(m,n);
gmax = maxc(amax);
/* loop through the tables: */
k = 0; cc = {};
do until k ge rows(t); k = k + 1;
   /* range of A cell: */ c = seqa(amin[k],1,amax[k]-amin[k]+1);
   /* nonzero multipliers in sum: */
      c = combin(n[k],c).*combin(t[k]-n[k],m[k]-c);
   if amin[k] gt 0; c = zeros(amin[k],1)|c; endif;
   if amax[k] lt gmax; c = c|zeros(gmax-amax[k],1); endif;
   cc = cc|c';
endo;
e = (cc .gt 0).*seqa(0,1,cols(cc))';
s = cc.*(r.^e);
retp(s./sumr(s));
endp;
