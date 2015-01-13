/* This proc computes basic descriptive statistics on a dataset
   with missing values supplied by user (not GAUSS codes).
   Missing values are counted and dropped from computations.
   Inputs: dat = data matrix if inp = 0,
                 file handle of closed gauss dataset if inp = 1
           inp = 0 if dat is data matrix,
                 1 if dat is file handle of closed gauss dataset
           mis = missing value code(s) --
                 single number if same for all columns of dat
                 column vector with values for each dat column otherwise
           names = names of columns of dat (0 if no printed output desired)
   Outputs: matrix with column means in row 1, sample SDs in row 2,
            mins in row 3, maxs in row 4, ranges in row 5,
            no. not missing in row 6, no. unique nonmissing in row 7.
*/
proc descmiss(dat,inp,mis,names);
local fname,dt,nc,nr,maxi,mini,ui,ret,k,nmis,cdat;
if inp;
   open dt = ^dat;
   dat = readr(dt,rowsf(dt));
   close(dt);
endif;
nc = cols(dat); nr = rows(dat);
if rows(mis) ne nc;
   if rows(mis) eq 1; mis = mis*ones(nc,1);
      else; "ERROR: No. missing values supplied to DESCMISS.G";
            "       doesn't match columns of data matrix."; end; 
   endif;
endif;
if 0$+names[1]; /* print summary table: */ print;
   if inp; "For dataset ";; $fname;;" - "; endif;
   "Summary of "$+ftocv(nr,1,0)$+" rows and "$+ftocv(nc,1,0)$+" columns:";
"              mean        SD       min        max       #obs #values   #miss";
endif;
k = 0; ret = zeros(nc,7);
do until k ge nc; k = k + 1; 
   cdat = delif(dat[.,k],dat[.,k] .eq mis[k]);
   ui = rows(unique(cdat,1)); nmis = nr-rows(cdat);
   if nmis lt nr;
      maxi = maxc(cdat); mini = minc(cdat); 
      ret[k,.] = meanc(cdat)~stdc(cdat)~mini~maxi~rows(cdat)~ui~nmis;
      else; ret[k,.] = -9*ones(1,6)~nmis;
   endif;
   if 0$+names[1]; /* print summary table: */
      format /rdn 2,0; $ftocv(k,1,0);;".";;
      format /lds 8,0; $names[k];; format /rds 9,2; ret[k,1:3];;
      format 10,2; ret[k,4];;  format 7,0; ret[k,5:7];
   endif;
endo; format 9,3;
retp(ret);
endp;
