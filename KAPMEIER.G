/* This proc computes the Kaplan-Meier survival probability for the times t
   without ordering the times.
   Inputs: t = failure or censoring times (individual data)
           d = failure indicator
           rep = repetiton count (0 if none)
   Outputs: f = unique failure times
            s = KM survival prob estimates at f
*/
proc (2) = kapmeier(t,d,rep);
local f,ns,nf,s;
if rows(rep) eq 1; rep = ones(rows(t),1); endif;
/* failure times: */ f = unique(selif(t,d),1);
/* no. at risk at f: */ ns = countge(t,f,rep);
/* no. failing at f: */ nf = counteq(t,f,d.*rep);
s = cumprodc((ns - nf)./ns);
retp(f,s);
endp;
