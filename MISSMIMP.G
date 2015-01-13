/* mean-imputes missing values in x */
proc missmimp(x,miscode);
local nc,mi,mu,j;
nc = cols(x);
miscode = ones(nc,1).*miscode;
mi = x .eq miscode';
mu = zeros(nc,1); j=0;
do until j ge nc; j = j+1;
   mu[j] = meanc(delif(x[.,j],mi[.,j])); 
endo;
retp(x-(x-mu').*mi);
endp;