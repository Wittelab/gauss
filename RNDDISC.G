/* This program generates discrete variates
   with distributions (1-sumr(p))~p over the range 0 to cols(p) */
/* Test program: */
let p[1,4] = .1 .3 .2 .25;
y = rnddisc(p[ones(10000,1),.]); format /rdn 10,3;
t = sumc((y .eq 0)~(y .eq 1)~(y .eq 2)~(y .eq 3)~(y .eq 4));
t'/10000; sumc(t);
end;
proc rnddisc(p);
retp(sumr(rndu(rows(p),1) .lt cumsumc(rev(p'))'));
endp;
