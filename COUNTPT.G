/* This proc counts the person time about each failure in failure-time data.
   Inputs: te = time of entry under observation
           tx = time of exit from observation
           d = status at exit = 1 if failure, 0 if right-censored
           a = vector of repetion counts (0 if none).
   Outputs: midp = unique failure times
            ucutp = upper points for intervals about midp
            lcutp = lower points for intervals about midp
            pt = person-time in intervals about midp.
*/
proc (4) = countpt(te,tx,d,a);
local n,nones,nin,ucutp,midp,lcutp,tintc,nin1,pt;
n = rows(tx);
/* Determine unique failure times and set up an interval about each: */
midp = unique(selif(tx,d),1); nin = rows(midp);
ucutp = (midp[2:nin] - midp[1:(nin-1)])/2; ucutp = ucutp|maxc(tx);
lcutp = minc(te)|ucutp[1:(nin-1)];
/* Count the person-time in each interval: */
nin1 = ones(nin,1); nones = ones(n,1); 
if rows(te) eq 1; te = te*ones(n,1); endif;
pt = minr((nones.*.ucutp)~(tx.*.nin1)) - maxr((nones.*.lcutp)~(te.*.nin1));
pt = maxr(pt~zeros(na,1));
if sumc(a);
   pt = sumc((a.*reshape(pt,n,nin))');
   else; pt = sumc(reshape(pt,n,nin)'); endif;
retp(midp,ucutp,lcutp,pt);
endp;