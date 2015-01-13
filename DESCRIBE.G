/* This proc computes basic descriptive statistics on a dataset
   that has no missing values.
   Inputs: dat = data matrix if inp = 0,
                 file name of closed gauss dataset if inp = 1
           inp = 0 if dat is data matrix,
                 1 if dat is file name of closed gauss dataset
           names = names of columns of dat (0 if no output desired)
   Outputs: matrix with column means in row 1, sample SDs in row 2,
            mins in row 3, maxs in row 4, ranges in row 5.
*/
proc describe(dat,inp,names);
local fname,dt,maxi,mini,ret,k;
if inp;
   let fname = "dat";
   open dt = ^fname;
   dat = readr(dt,rowsf(dt));
   close(dt);
endif;
maxi = maxc(dat); mini = minc(dat);
ret = meanc(dat)~stdc(dat)~mini~maxi~(maxi-mini);
if 0$+names[1]; /* print summary table: */ print;
   if inp; "For dataset ";; $fname;;" - "; endif;
   "Summary of ";; format 7,0;
    rows(dat);;" rows and ";; format 3,0; cols(dat);; " columns:";
"                     mean     std dev         min         max       range     #vals";
k = 0;
do until k ge rows(ret); k = k + 1;
   format /rdn 3,0; $ftocv(k,1,0);;". ";;
   format /lds 8,0; $names[k];; format /rds 11,3; ret[k,.];; 
   format 9,0; rows(unique(dat[.,k],1));
endo; format 9,3;
endif;
retp(ret);
endp;
