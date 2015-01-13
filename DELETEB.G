/* This proc approximates the delta-betas (observation-deleted betas
   MINUS fitted betas) for models linear in the natural parameter of a
   one-parameter exponential-family error distribution.
   for generalized-linear (in the natural parameter) regressions.
   Inputs: x = regressor (design) matrix, including constant if present
           e = observed minus fitted y (residual on mean scale)
           w = final fitted weight
           b = coefficients
           bcov = estimated covariance matrix for b
   Output: matrix of approximations to one-step delta-betas
           -- row i is approx. delta-beta for deleting observation i.
 */
proc deleteb(x,e,w,b,bcov);
local crwx,ch;
    /* h = sqrt(w).*x*bcov*(x.*sqrt(w))' may be too big to fit in memory, so
       use the following trick to avoid computing h. h = crwx'crwx but we
       don't need all of h, just its diagonal, which is sumc(crwx.*crwx): */
       crwx = chol(bcov)*(sqrt(w).*x)';
    /* variance adjustment factor for approximate delta-betas
       (needed because bcov is not recomputed for each observation): */
        ch = 1-sumc(crwx.*crwx);
/* return approx. delta-betas (each row is an approximation to the
   first Newton step toward what b would become after dropping
   the corresponding observation): */
retp(-(e./ch).*(x*bcov));
endp;
