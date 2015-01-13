/* This proc returns the generalized logistic transform of x given s,
   a matrix of nonzero proportions each row of which sums to less than 1
   (useful for multinomial logistic modelling).
   NOTE: x and s' must be conformable for x.*s' */
proc expitg(x,s);
local e;
e = exp(x.*s');
retp(e./(1+sumr(e)));
endp;
