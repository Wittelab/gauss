/* Returns mean of a normal(0,1) random variable Z given a<Z<b, computed by quadrature.
   Inputs:  a = col vector of lower limits
            b = col vector of upper limits
   Outputs: col vector of E(Z|a<Z<b)
   Globals: _intord = order of accuracy of quadrature 
           (default = 40,the maximum provided by Gauss; values below 12 are severely inaccurate in test examples)
 */
/* Test program: */ output file = NORMX.OUT reset; format 12,4;
p = .025|.05|seqa(.1,.1,9)|.95; z = invnorm(1-p);
"For given power: Z cutoff for that power,";
" bias factors in normal estimators selected by 2-sided .05 significance (tsig)";
" upper .025 significance (usig), 1/bias for 2-sided .05 nonsignificance(nsig)";
"upper power, 2-tail power,         z,    bias tsig,    bias usig,  1/bias nsig -"; /* mad sig, mad nonsig not implemented */
bn = normx(z-3.92,z); bl = normx(-6,z-3.92); bu = normx(z,6);
/* for bias tsig: */ wl = 1-cdfnc(z-3.92); wu = cdfnc(z); pow = wl+wu; bt = (wl.*bl+wu.*bu)./pow;
"for variance = ";; v = 1/8; v;; ":";;
p~pow~z~exp((bt~bu~(-bn))*sqrt(v));
end; 
proc normx(a,b);
local r;
r = ones(1,maxc(rows(a)|rows(b))); a = a'.*r; b = b'.*r; _intord = 40;
retp(intquad1(&znd,b|a)./(cdfnc(a)-cdfnc(b))');
endp;
/* z*standard normal density: */
proc znd(z); 
retp( z./sqrt(2*pi*exp(z.*z)) );
endp;