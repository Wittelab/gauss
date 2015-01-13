/* This proc returns indices of modes (local maxima) in a vector f;
   returns empty vector {} if none: */
proc modeindc(f);
local n,mi,i; n = rows(f);
   if f[1] gt f[2]; mi = 1; else; mi = {}; endif;
   i = 1;
   do until i ge n-1; i = i + 1;
      if f[i] gt f[i-1] and f[i] gt f[i+1]; mi = mi|i; endif;
   endo;
   if f[n] gt f[n-1]; mi = mi|n; endif;
retp(mi);
endp;
