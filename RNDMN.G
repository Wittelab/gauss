/* Test program (should print back a number near r): */
np = 3; m = zeros(1,np); r = .6;
r = corrx(rndmn(m,diagrv(r*ones(np,np),ones(np,1)),25000));
(sumall(r)-np)/(np*(np-1));
end;
/* This proc generates a matrix with multivariate normal rows.
   Inputs: m = mean row vector for variates
           c = covariance matrix for variates
           n = number of output rows (no. of variates)
   Output: x = matrix of multinormal(m,c) variates
*/
proc rndmn(m,c,n);
retp(m + rndn(n,cols(m))*chol(c));
endp;
