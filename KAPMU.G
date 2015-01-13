/* This proc computes the Kaplan-Meier survival probability for the times t
   without ordering the times.
   Inputs: t = unique failure times
           a = no. of failures at t (may be set to 0 if always 1 failure)
   Outputs: s = KM survival prob estimates
 */
proc kapmu(t,a);
local n,ns,lp,s,i;
n = rows(t);
ns = countge(t,t,a);
if sumc(a);
   lp = ln((ns - a)./ns);
   s = ones(n,1);
   do until i ge n; i = i + 1;
       s[i] = exp(sumc((t[i] .ge t).*lp));
   endo;
   else;
       s = ((ns-1)./ns).^(n - ns + 1);
endif;
retp(s);
endp;
