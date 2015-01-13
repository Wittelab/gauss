/* Test program:
t0 = date;
i = 0; t = 400; p = 0; r = 400; rt = r*t;
do until i ge t; i = i + 1;
   p = p + rows(unique(rndsrep(r,r),1))/rt;
endo;
"Total run time: ";; etstr(ethsec(t0,date));
format /rdn 7,4; "Average proportion selected:";; p;
"(stochastically converges to ";; 1-exp(-1);; " as t and r increase)";
end;
*/
/* This proc selects a random sample of size s from the numbers 1 to n,
   with replacement.
   NOTE: if r = rows(x), x[rndsrep(r,r),.];
         produces a bootstrap sample from the data matrix x.
*/
proc rndsrep(n,s);
retp(ceil(n*rndu(s,1)));
endp;
