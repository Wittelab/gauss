/* This proc computes the WLS coefficients and 
   inverse sum-of-covariate-cross-products matrix from x,y,w.
 Inputs: x = design matrix (MUST INCLUDE CONSTANT IF ONE IS DESIRED),
         y = outcome vector,
         w = weight vector (set to 1 for equal-weight OLS regression).
 Outputs: b = coefficient estimates,
          iscp = inverse of sum-of-covariate-cross-products matrix
                 (this is the inverse-information covariance matrix for b
                  if w is the inverse-variance vector for y).
*/
proc (2) = wlsquick(x,y,w);
local xw,b,iscp;
xw = (w.*x)';
iscp = invpd(xw*x);
b = iscp*xw*y;
retp(b,iscp);
endp;
