/* This procedure computes 95% limits, medians, and two-sided p values
   from a Monte-Carlo distribution.
   Inputs: b = distribution matrix (rows are observations)
           names = names of columns of b (0 if no printed output desired)
   Outputs: lim = lower limit | upper limit
            med = medians
            p = two-sided p value
*/
proc (3) = mclimits(b,names);
local nr,p25,p975,k,lim,med,bsort,p,up;
nr = rows(b);
p25 = maxc(floor(.025*nr)|1); p975 = minc(ceil(.975*nr)|(nr-p25));
if 0$+names[1];
   " 2.5th, 50th, and 97.5th percentiles of Monte Carlo distributions,";;
   " & p values:";
   if nr%40 ne 0;
      " WARNING: number of trials not divisible by 40 -";
      "  lower & upper CL percentiles are ";;
      format 6,4; p25/nr;;" and ";; p975/nr;
   endif;
endif;
k = 0; lim = {}; med = {};p = {};
do until k ge cols(b); k = k + 1;
   bsort = sortc(b[.,k],1);
   lim = lim|bsort[p25 p975]';
   med = med|median(bsort);
   up = sumc(bsort .ge 0)/nr;
   p = p|(2*minc(up|1-up));
   if 0$+names[1];
      format /ldn 10,0; $names[k];;
      format /rdn 10,4; lim[k,1]~med[k]~lim[k,2];;",";;
      format 12,4; p[k]; endif;
endo;
retp(lim,med,p);
endp;
