/* This proc computes some basic descriptive statistics for the
   distributions of the columns of x, including percentiles useful for
   analysis of Monte-Carlo simulation results. Will also supply histograms
   to accompany these statistics. See DESCRIBE.G for a simpler, faster
   running proc to describe actual data, including GAUSS datasets.
   Inputs: x = data matrix
           histo = number of bins for histograms (set to 0 for no histograms)
           names = names for column of x (set to 0 if no printed ouput desired)
   Outputs: matrix for which row k contains mean, standard deviation, and
            median, and upper & lower 2.5th percentiles of x[.,k] (these
            percentiles are exact only if rows(x) is divisible by 40).
*/
proc desc(x,histo,names);
local nr,lp25,up25,est,k,ord,ret;
if histo; library pgraph; endif;
nr = rows(x);
/* multiplier for standard deviation: */ let est = 0 -1.96 1.96;
/* rank of 2.5th percentile: */ lp25 = 1 + floor(nr/40);
/* rank of 97.5th percentile: */ up25 = ceil(nr*39/40);
if names[1]; /* print table header: */
   format /rdn 7,0;"Summary of ";; nr;; " rows and ";;
   format 3,0; cols(x);; " columns";
" (if no. rows is not divisible by 40, lower 2.5th %tile is largest value with";
"  at least 97.5% of values at or above it, and upper 2.5th %tile is smallest";
"  value with at least 97.5% of values at or below it):";
"              m=mean  m-1.96*SD  m+1.96*SD;     median  lower 2.5% upper 2.5%";
endif;
k = 0; ret = {};
do until k ge cols(x); k = k + 1;
   /* sort in ascending order to find percentiles: */ ord = sortc(x[.,k],1);
   ret = ret|meanc(x[.,k])~stdc(x[.,k])~median(x[.,k])~ord[lp25 up25]';
   if 0$+names[1]; /* print row k of table: */
      format /ldn 9,0; $names[k];; format /rdn 11,4;
      (ret[k,1] + est'ret[k,2]);;";";;ret[k,3]~ord[lp25 up25]';
   endif;
   if histo; /* give histograms: */
      output off; "Press enter to see histogram of ";; $names[k];
      wait; call hist(x[.,k],histo); output on;
   endif;
endo;
retp(ret);
endp;
