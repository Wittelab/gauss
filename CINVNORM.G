/* This procedure finds the value that cuts off the upper p proportion of a
   standard normal distribution F(x), i.e., solves p = 1 - F(x) for x, using
   the Derenzo approximation (see Hoaglin, Am Statist 1989; 43:289).
   Absolute error less than .00013 for p > 10E-7.
   NOTE THAT GAUSS 3.2.13 & LATER HAS A EQUALLY ACCURATE PROCEDURE
   1-CDFNI(x). */
proc cinvnorm(p);
local d,y;
d = (p .gt .5);
y = -ln(2*(p - d.*(2*p - 1)));
y = y.*y.*(y.*(4*y + 100) + 205)./(y.*(y.*(2*y + 56) + 192) + 131);
y = sqrt(y).*((-1).^d);
retp(y);
endp;

