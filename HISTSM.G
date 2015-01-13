/* This proc computes the smoothed histogram and histogram modes of a vector x
   using an Epanechnikov kernel (rec. for smooth underlying densities only).
   See DENSITY.G for a simpler density-estimation proc.
   Inputs: x = data vector
           rep = repetition counts (0 if none)
           h = radius of window of nonzero weights (smoothing parameter)
               -- 1/2 standard deviation of x can be a good starting value
           n = no. of histogram points
           name = name of x (0 for no printed output)
           dgraph = 1 for graph of histogram estimate f
   Outputs: f = smoothed percentile histogram using n evenly spaced points
            m = vector of estimated modes
                (if too many, increase h &/or decrease n)
   WARNING: A mode of f that extends across 2 or more values of f
            (i.e., a mode occurring at a flat part of f) will be missed!
*/
proc (2) = histsm(x,rep,h,n,name,dgraph);
local l,u,d,e,f,ifm,fm,m,mi,i,ni;
if sumc(rep) eq 0; rep = 1; endif;
l = minc(x); u = maxc(x); d = (u-l)/(2*(n-1));
if d eq 0; "HISTSM.G: NO VARIATION IN INPUT."; retp(-9999,-9999); endif;
/* evaluation points: */ e = seqa(l,2*d,n);
/* smoothed histogram: */
   f = sumc(rep.*max(1 - ((e' - x)./h)^2,0)); f = 100*f/sumc(f);
   ifm = maxindc(f); fm = f[ifm]; mi = modeindc(f);
if dgraph; library pgraph; graphset; _pdate =""; fonts("microb");
   _plctrl = 0; _plwidth = 2; _pltype = 6; ytics(0,fm,ceil(100*fm)/100,4);
   xlabel(name); ylabel("Smoothed histogram"); xy(e,f);
endif;
if 0$+name; format /rds 6,4;
   "For "$+name$+" using "$+ftocv(n,1,0)$+" points from ";;
    l;; "to " u;;" with kernel radius ";; h;;":"; format 8,4;
   " Histogram modes at         :";; e[mi]';
   " Histogram height at modes  :";; f[mi]';
   " Histogram maximum at       :";; e[ifm];
   " Histogram height at maximum:";; fm;
   " Histogram minimum at       :";; e[minindc(f)];
   " Histogram height at minimum:";; minc(f);
   " Numeric integral           :";; ni = 2*d*sumc(f); ni;; 
   if ni lt .99 or ni gt 1.01; "<=WARNING, should be 1"; else; print; endif;
endif;
retp(f,e[mi]);
endp;
