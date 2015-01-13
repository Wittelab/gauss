/* This proc computes the Epanechnikov kernel density estimate for x.
   Suitable only for smooth underlying densities.
   See HISTSM.G for a more detailed proc that graphs f and finds its modes.
   Inputs: x = input vector
           h = radius of window of nonzero weights
           xp = SORTED points on x axis at which to estimate the density
   Output: f = density estimate for x at the points xp.
*/
proc density(x,h,xp);
local f;
f = rows(x) - sumc(min(((xp' - x)./h)^2,1));
retp(f./sumc(f));
endp;
