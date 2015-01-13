/* This procedure generates an m by n matrix of 
   trapezoidal variates with changepoints a,b,c,d.
   For triangular(m,n,a,c,d), set b=c.
*/
/* test program -- note edge effects in graphs:
t = rndtrapz(5000,1,1,2,3,5);
call histsm(selif(t, t .lt 2),0,0.05,500,"trapz",1);
call histsm(selif(t, (t .ge 2).*(t .le 3)),0,0.05,500,"trapz",1);
call histsm(selif(t, t .gt 3),0,0.05,500,"trapz",1);
end; */
proc rndtrapz(m,n,a,b,c,d);
local t,u,v;
t = (a + b + (c+d-a-b).*rndu(m,n))/2;
t = (t .lt b).*(a + sqrt((b-a).*(2*t-a-b)))
    + (t .ge b).*(t .le c).*t
    + (t .gt c).*(d - sqrt(2*(d-c)*(t-c)));
retp(t);
endp;
