/* This proc computes the Mantel-Haenszel trend statistic for 2xJxK arrays.
   Inputs: a = case-count array (row 1 of each stratum), one ROW per stratum
           n = case+noncase count array (denominators), one ROW per stratum
           s = column vector of scores for columns of a and n.
   Output: Mantel-Haenszel chi statistic for trend (NOT CHI-SQUARED).
*/
proc mhtrend(a,n,s);
local sa,sn,sb,as,es,v0,z;
sa = sumr(a); sn = sumr(n); sb = sn-sa;
as = a*s; es = (n*s)./sn;
v0 = (bn*(bs^2)./sn - es^2)'(sa.*sb./(sn-1));
retp((as - sa'es)/sqrt(v0));
endp;
