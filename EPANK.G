/* This proc computes the Epanechnikov kernel matrix for a vector of values
   using a fixed bandwidth (see Silverman, 1986).
   Inputs: x = input vector (usually a vector of regressor values)
           h = bandwidth (radius for positive weights)
   Output:  k = kernel weight matrix -- each row contains the weights
             assigned to each element of x for the kernel centered at
             that element of x. WEIGHTS ARE NORMALIZED TO SUM TO 1.
   WARNING: k can be a huge matrix - it has rows(x) rows and columns!
*/
proc epank(x,h);
local k;
k = max(1 - ((x - x')./h)^2,0);
k = k./sumc(k');
retp(k);
endp;