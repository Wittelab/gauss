/* This proc does sensitivity analysis of the impact of misclassification
   of columns (within rows) on the Mantel-Haenszel odds ratio for
   K 2x2 tables, allowing for possible differential errors across rows.
   Inputs: a1,a0 = first-row entries of observed 2x2 tables
           b1,b0 = second-row entries of observed 2x2 tables
           se,sp = H by 2 matrices of sensitivities and specificities for the
                   first then second rows of the tables (these are assumed
                   constant across strata); H corrected odds ratios
                   will be generated, one for each row of se,sp.
           rname = name of row variable (if 0, no printed output)
           cname = name of column variable (if 0, no printed output)
           snames = names of stratification variables (ignored if 0)
   Outputs: orc = vector of MH odds ratios computed from the corrected counts,
                  one for each row of se,sp
*/
proc (3) = misclass(a1,a0,b1,b0,se,sp,rname,cname,snames);
local h,fp,detr,i,orc,a1c,b1c,lse,p,headc;
h = rows(se);
/* false-positive probs, & determinants: */
   fp = 1 - sp; detr = se.*sp - (1-se).*fp;
/* will be corrected odds-ratio vector: */ orc = -ones(h,1); i = 0;
do until i ge h; i = i + 1;
   if prodr(detr[i,.] gt 0);
   /* corrected counts: */ a1c = (sp[i,1]*a1 - fp[i,1]*a0)/detr[i,1];
      b1c = (sp[i,2]*b1 - fp[i,2]*b0)/detr[i,2];
      { orc[i],lse,p } = ormh2x2(a1c,a1+a0-a1c,b1c,b1+b0-b1c,0);
   endif;
endo;
if rname[1] and cname[1]; print;
   " PROC MISCLASS.G: Sensitivity of the odds ratio to misclassification.";
   "      Column variable is "$+cname$+", row variable is "$+rname$+".";
   if snames[1];
      " Corrected Mantel-Haenszel odds ratio with stratification by";
      format 8,0; $snames';;".";
   endif;
   "      Se = sensitivity = true-positive "$+cname$+" probability.";
   "      Sp = specificity = true-negative "$+cname$+" probability.";
   format 20,2; let headc = "ROW 1  "  "ROW 2  "; $headc';;
   "   Corrected OR";
   format 10,2; let headc = "Se"  "Sp"; $headc';;$headc';
   se[.,1]~sp[.,1]~se[.,2]~sp[.,2]~orc;
   " A corrected OR of -1 indicates the classification probabilities in";
   " that row yield negative cell counts and so are not compatible with";
   " the observed counts.";
endif;
retp(orc,a1c,b1c);
endp;
