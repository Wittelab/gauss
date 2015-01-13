/* Test: invnorm(.05|.95); */
/* This procedure finds the value that cuts off the lower p proportion of a
   standard normal distribution F(x), i.e., solves p = F(x) for x, using
   the Derenzo approximation (see Hoaglin, Am Statist 1989; 43:289).
   Absolute error less than .00013 for p > 10E-7.
   NOTE THAT GAUSS 3.2.13 & LATER HAS A EQUALLY ACCURATE PROCEDURE
   CALLED CDFNI(x). */
proc invnorm(p);
local y;
y = -ln(2*(p - (p .gt .5).*(2*p - 1)));
y = y.*y.*(y.*(4*y + 100) + 205)./(y.*(y.*(2*y + 56) + 192) + 131);
y = -sqrt(y).*((-1).^(p .gt .5));
retp(y);
endp;

