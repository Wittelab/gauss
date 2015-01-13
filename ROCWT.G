/* This proc does weighted ROC (receiver operating characteristic) curve analysis.
   Inputs: y = true binary outcome status
           py = fitted or predicted values for y
           wt = weights for each observation (e.g., delta/k, where delta is an
                indicator of whether desired exposure level used and k is probability 
                given covariates of desired exposure level)
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
proc (5) = rocwt(y,py,wt,yname,groc);
/* number of intervals */ declare matrix _ni = 40;
local a,nc,cuts,pvp,pvn,tp,fp,i,pos,npos,nneg,area;
/* no. positives and no. cutpoints: */ a = sumc(wt.*y); nc = 1 + _ni;
/* cutpoints: */ cuts = seqa(0,1/_ni,nc);
/* positivity indicators and no. positive using each cutpoint: */
   pos = py .gt cuts'; npos = sumc(pos.*wt); nneg = sumc((1-pos).*wt);
/* true & false positive proportions: */
   tp = (pos'(y.*wt))/a; fp = (pos'((1-y).*wt))/(sumc(wt.*(1-y)));
/* positive & negative predictive values: */
   pvp = (pos./miss(npos',0))'(y.*wt);
   pvn = ((1-pos)./miss(nneg',0))'((1-y).*wt);
   area = (firstdif(fp)'(tp[1:_ni]+tp[2:nc]))/2;
if (0$+yname); /* print results: */ format /rdn 11,3;
   "  Cutpoints and corresponding true & false positive proportions and";
   "  positive and negative predictive values, as computed from the";
   "  fitted response probabilities (a dot means undefined):";;
   cuts~tp~(1-tp)~fp~(1-fp)~pvp~pvn;
   print; "  Area under the fitted ROC curve";; area;
   if groc; library pgraph; title("ROC curve for "$+yname);
      xlabel("False-positive proportion (1-specificity)");
      ylabel("True-positive propotion (sensitivity)"); xy(fp,fp~tp);
   endif;
endif;
retp(pvp,pvn,tp,fp,area);
endp;


