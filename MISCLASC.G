/* This proc does misclassification correction
   within rows (a,b; c,d)) on the crude odds ratio for
   2x2 tables, allowing for possible differential errors across rows.
   Inputs: a,b = first-row entries of observed 2x2 tables
           c,d = second-row entries of observed 2x2 tables
           se,sp = H by 2 matrices of sensitivities and specificities for the
                   first then second rows of the tables (these are assumed
                   constant across strata); H corrected odds ratios
                   will be generated, one for each row of se,sp.
           rname = name of row variable (if 0, no printed output)
           cname = name of column variable (if 0, no printed output)
   Outputs: orc = vector of crude odds ratios computed from the corrected counts,
                  one for each row of se,sp
            ac,bc,cc,dc = corrected counts
*/
proc (5) = misclass(a,b,c,d,se,sp,rname,cname);
local fp,orc,ac,bc,cc,dc,headc;
/* false-positive probs: */ fp = 1-sp;
   /* corrected counts and OR: */ 
      ac = (a - fp[.,1].*(a+b))./(se[.,1]-fp[.,1]); bc = a+b-ac;
      cc = (c - fp[.,2].*(c+d))./(se[.,2]-fp[.,2]); dc = c+d-cc;
   orc = ac.*dc./(bc.*cc);
if rname[1] and cname[1]; print;
   " PROC MISCLASC.G: misclassification corrction.";
   "      Column variable is "$+cname$+", row variable is "$+rname$+".";
   "      Se = sensitivity = true-positive "$+cname$+" probability.";
   "      Sp = specificity = true-negative "$+cname$+" probability.";
   format 20,2; let headc = "ROW 1  "  "ROW 2  "; $headc';;
   "   Corrected OR";
   format 10,2; let headc = "Se"  "Sp"; $headc';;$headc';
   se[.,1]~sp[.,1]~se[.,2]~sp[.,2]~orc;
   " A negative corrected OR indicates the classification probabilities in";
   " that row yield negative cell counts and so are not compatible with";
   " the observed counts.";
endif;
retp(orc,ac,bc,cc,dc);
endp;
