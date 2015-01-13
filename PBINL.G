/* Test program:  format 9,4; t = 10000; n = 4*ones(3*t,1);
   let k = 0 2 4; k = k.*.ones(t,1);
   let p = 0 .5 1; p = ones(t,1).*.p;
t0 = date; call pbinl(k,p,n);
"Total run time pbinl: ";; etstr(ethsec(t0,date));
t0 = date; call pbinf(k,p,n);
"Total run time pbinf: ";; etstr(ethsec(t0,date));
/*  k~p~n~pbinl(k,p,n)~pbinf(k,p,n); */
 end; */
/* This proc computes the exact probability of a binomial(p,n) outcome k
   using log factorials. Tends to be faster than pbinf.g if 0<p<1.
   Inputs: k = Binomial outcome (integer between 0 and n);
           p = Bernoulli probability
           n = no. of Bernoulli trials
   Output: Pr{bin(p,n) .eq k} */
proc pbinl(k,p,n);
local test0,test1;
test0 = p .eq 0; test1 = p .eq 1;
retp(test0.*(k .eq 0) + test1.*(k .eq n) + (1-test0).*(1-test1).*
     exp(lnfact(n)-lnfact(k)-lnfact(n-k)+k.*ln(p+test0)+(n-k).*ln(1-p+test1)));
endp;
