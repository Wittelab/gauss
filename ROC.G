/* test program:  let y = 1 1 1 1 1 0 0 0 0 0;
let py = .7 .3 .5 .5 .4 .5 .4 .6 .2 .1;
_ni = 10; call roc(y,py,"y",1); end;
load dat[130,8] = \gauss\prg\testdat\rosm3.dat;
y = dat[.,1:4]; py = dat[.,8 5:7]; let names = g c2 c5 c11;
i = 2; _ni = 60; call roc(y[.,i],py[.,1],names[i],0); end; */
/* This proc does ROC (receiver operating characteristic) curve analysis.
   Inputs: y = true binary outcome status
           py = fitted or predicted values for y
           yname = name of y (set to 0 for no printed output or graph)
           groc = 1 for graph of ROC curve, 0 for none (requires yname)
                  -- fineness of graph (smoothness) is controlled by
                     the global _ni described below
   Outputs: pvp = predicted value positive at each graph point
            pvn = predicted value negative at each graph point
            tp = true positive proportion (sensitivity) at each point
            fp = false positive proportion (1 - specificity) at each point
            area = area under roc curve
   Global: _ni = number of intervals between plotted cutpoints (default 40)
                 -- set low (under 20) for more smooth plot,
                    set high (over 40) for finer detail in plot.
*/
proc (5) = roc(y,py,yname,groc);
/* number of intervals */ declare matrix _ni = 40;
local a,nc,cuts,pvp,pvn,tp,fp,i,pos,npos,area;
/* no. positives and no. cutpoints: */ a = sumc(y); nc = 1 + _ni;
/* cutpoints: */ cuts = seqa(0,1/_ni,nc);
/* positivity indicators and no. positive using each cutpoint: */
   pos = py .gt cuts'; npos = sumc(pos);
/* true & false positive proportions: */
   tp = (pos'y)/a; fp = (pos'(1-y))/(rows(y)-a);
/* positive & negative predictive values: */
   pvp = (pos./miss(npos',0))'y;
   pvn = ((1-pos)./miss(rows(y)-npos',0))'(1-y);
if (0$+yname); /* print results: */ format /rdn 11,3;
   "  Cutpoints and corresponding true & false positive proportions and";
   "  positive and negative predictive values, as computed from the";
   "  fitted response probabilities (a dot means undefined):";;
   cuts~tp~(1-tp)~fp~(1-fp)~pvp~pvn;
   area = (firstdif(fp)'(tp[1:_ni]+tp[2:nc]))/2;
   print; "  Area under the fitted ROC curve";; area;
   if groc; library pgraph; title("ROC curve for "$+yname);
      xlabel("False-positive proportion (1-specificity)");
      ylabel("True-positive propotion (sensitivity)"); xy(fp,fp~tp);
   endif;
endif;
retp(pvp,pvn,tp,fp,area);
endp;


