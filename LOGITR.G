/* This proc calculates risks, risk ratios, and risk differences and
   their confidence limits from logistic coefficients and their covariances,
   for the indicated regressor vectors.
   Inputs: b = coefficients
           cov = covariance matrix of coefficients
           e = vector of index regressor values
           r = vector of reference regressor values
           const = 1 if constant in b, 0 otherwise
           xnames = names of regressors
   Outputs: pe = index logistic-risk estimate & 95% ci
            pr = ref logistic-risk estimate & 95% ci
            rr = risk-ratio estimate & 95% ci
            rd = risk-difference estimate & 95% ci
*/
proc (4) = logitr(b,cov,em,rm,const,xnames);
declare matrix _z = 1.96;
local n,i,e,r,er,lee,lr,pe,pr,ve,vr,erad,rrad,d,cer,lrr,rlrr,rd,rrd,xi;
n = cols(em);
if const; em = ones(1,n)|em; rm = ones(1,n)|rm; endif;
i = 0;
do until i ge n; i = i + 1;
   e = em[.,i]; r = rm[.,i]; er = e~r;
   lee= b'e; lr = b'r;
   ve = e'cov*e; vr = r'cov*r;
   erad = _z*sqrt(ve); rrad = _z*sqrt(vr);
   pe = expit(lee|(lee-erad)|(lee+erad));
   pr = expit(lr|(lr-rrad)|(lr+rrad));
   /* derivative of logistic transforms: */
   d = diagrv(eye(2),(pe[1]*(1-pe[1]))|(pr[1]*(1-pr[1])));
   /* covariance matrix of pe,pr: */ cer = d'er'cov*er*d;
   /* variance of pe, pr: */ ve = cer[1,1]; vr = cer[2,2];
   /* derivative of risk-ratio transform: */
      d = (1/pr[1])|(-pe[1]/(pr[1]*pr[1]));
   /* log risk-ratio: */ lrr = ln(pe[1]/pr[1]);
   /* log risk-ratio radius: */ rlrr = _z*sqrt(d'cer*d)/(pe[1]/pr[1]);
   rr = exp(lrr|(lrr-rlrr)|(lrr+rlrr));
   /* risk difference: */ rd = pe[1] - pr[1];
   /* risk-difference radius: */ rrd = _z*sqrt(ve + vr - 2*cer[1,2]);
   rd = rd|(rd-rrd)|(rd+rrd);
   xi = seqa(const+1,1,rows(e)-1);
   "For regressors ";; $xnames';
   "at index level ";; e[xi]';
   "and ref level  ";; r[xi]';
   " - point estimates and confidence limits for index risk,";
   "   reference risk, risk ratio, & risk difference:";
    format /rdn 15,5; pe'; pr'; rr'; rd'; print;
endo;
retp(pe,pr,rr,rd);
endp;

