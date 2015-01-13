/* This procedure checks all interactions among covariates with indices
   indint.
   Inputs: x = design matrix
           e = scores (observed - expecteds) for reduced model
           w = weights for reduced model
           indint = indices of x columns to test for interaction
           names = names of x columns to test for interaction
   Outputs: sstata = score statistic for all interactions
            sstate = vector of score statistics for each interaction
 */
proc (2) = scoreint(x,e,w,indint,names);
local xint,sstata,sstate,pval,nameint,k;
 "Score statistics and p values for interactions:"; format 10,4;
          xint = inter(x[.,indint]);
          "  Test of all specified interactions:";
          { sstata,pval } = scoretst(x~xint,cols(xint),e,w,"all");
          "  Tests of each specified interaction:";
          /* Name the interactions: */ nameint = ninter(names);
          sstate = zeros(cols(xint),1); k = 0;
          do until k ge cols(xint); k = k + 1;
             { sstate[k],pval } = scoretst(x~xint[.,k],1,e,w,nameint[k]);
          endo;
retp(sstata,sstate);
endp;
