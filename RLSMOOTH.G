/* This proc computes the smoothed running-weighted linear regression of y on x.
   Inputs: x = regressor
           y = regressand
           w = weight for y (if scalar, no weights will be used)
               -- if y is a proportion or rate, w should be the denominator
           hp = half-bandwidth (radius for positive weights, as proportion of
               range of x)
           &kern = &f, where f is the user-specified kernel function
           splot = 1 if plot of smooth is desired, 0 otherwise
           xname = name of x (0 for no output if plot = 0)
           yname = name of y (0 for no output if plot = 0)
           sboot = number of bootstrap plots of smooth desired
    Output: s = smoothed (predicted) values of y
*/
proc  rlsmooth(x,y,w,hp,&kern,splot,xname,yname,sboot);
local kern:proc,nr,h,k,mx,my,dx,b,s,ind,e,sm,indices,i,iboot,pct,sb;
nr = rows(y); h = hp*rangec(x);
if rows(w) gt 1; k = w'.*kern(x,h); k = k./sumc(k');
   else; w = ones(nr,1); k = kern(x,h); endif;
/* running-weighted means: */ mx = k*x; my = k*y;
/* running-weighted x-deviations: */ dx = x - mx;
b = k*(dx.*dx);
if sumc(b .eq 0);
   "zero divide - cannot smooth this sample using given bandwidth.";
   retp(0);
endif;
/* running-weighted coefficients: */ b = (k*(dx.*(y - my)))./b;
/* running-weighted smooth: */ s = my + b.*dx;
/* recenter the smooth: */ s = s + (w/sumc(w))'(y - s);
/* residual: */ e = y - s;
if xname[1]; print;
   "PROC RLSMOOTH.G: RUNNING WEIGHTED REGRESSION-LINE SMOOTHER";
   "                 regressing "$+yname$+" on "$+xname$+".";
   "Number of points to smooth         : ";; format 6,0; nr; format 9,3;
   "Weighted residual sum (should be 0):";; w'e;
   "Weighted residual sum of squares   :";; w'(e.*e);
   "Trace of smoother matrix           :";; tracemat(k); format 8,1;
   "Bandwidth as percent of range      :";; 200*hp;
endif;
if splot or sboot; library pgraph; xlabel(xname); ylabel(yname);
   ind = sortind(x); endif;
if splot; /* plot the smoother: */ title("Running weighted line smooth");
   let _plctrl = 1 -1; let _psymsiz = 1 2; let _pstype = 8 14; _pltype = 6;
   xy(x[ind],s[ind]~y[ind]);
endif;
if sboot; /* bootstrap the residuals & refit the smoother: */ clear k;
   sm = s[ind]~zeros(nr,sboot);
   indices = seqa(1,1,nr); i = 0;
   do until i ge sboot; i = i+1;
      iboot = ceil(nr*rndu(nr,1));
      pct = 100*rows(unique(iboot,1))/nr;
      format /rdn 4,0; "Percent selected for bootstrap sample ";; i;;
      ":";; format 5,1; pct;
      sb = rlsmooth(x,s+e[iboot],w,hp,&kern,0,0,0,0);
      if rows(sb) eq 1; i = i-1; continue; endif;
      sm[.,i+1] = sb[ind];
   endo;
   format /rdn 5,0;
   title("Smooth & "$+ftocv(sboot,1,0)$+" bootstrap replicates");
   _psymsiz = 1|(4*ones(sboot,1)); _plctrl = 1|zeros(sboot,1);
   _pltype = 6|(2*ones(sboot,1)); xy(x[ind],sm);
endif;
retp(s);
endp;
