/* This proc produces a matrix, a, which converts a coefficient vector beta
   wtih regression-coded indicators (j-1 variables coded 0,1 for j categories)
   to a coefficent vector alpha with ANOVA dummy variables (j variables coded
   -1,1 for j categories); also supplied is a matrix, b, which converts alpha
   back to beta.
   a and b are conditional inverses of one another, i.e.,
           a*b*a = a and b*a*b = b,
   but are not Moore-Penrose inverses of one another (a*b is not symmetric).
   NOTES: 1. Indicator groups MUST be contiguous in beta!
          2. Indicator groups MUST form the terminal segment of beta!
          3. beta and alpha MUST have intercept in first place!
   Inputs: s = vector of indices of starts of indicator groups in beta
           k = rows(beta)
   Outputs: a = regression-to-ANOVA matrix
            b = ANOVA-to-regression matrix
            t = vector of indices of starts of indicator groups in b
*/
proc (3) = reganova(s,k);
local n,j,j1,i,c,d,a,b,za,zb,t;
n = rows(s);
j = (((s[2:n]-1)|k) - s + 1);
t = s + seqa(0,1,n);
a = zeros(j[1],s[1]-1);
b = zeros(j[1]-1,s[1]-1);
c = {}; d = {}; za = {}; zb = {};
i = 0;
do until i ge n; i = i + 1;
   j1 = j[i] - 1;
   c = c~(ones(1,j1)/j[i]);
   d = d~1~(-ones(1,j1));
   if i gt 1;
      za = zeros(rows(a),j1);
      zb = zeros(rows(b),j[i]);
      a = a|zeros(j[i],cols(a));
      b = b|zeros(j1,cols(b));
   endif;
   a = a~(za|((-ones(1,j1)/(2*j[i]))|(eye(j1)/2-ones(j1,j1)/(2*j[i]))));
   b = b~(zb|(2*((-ones(j1,1)~eye(j1)))));
endo;
if s[1] gt 2;
   a = (eye(s[1]-1)~(c|zeros(s[1]-2,cols(c))))|a;
   b = (eye(s[1]-1)~(d|zeros(s[1]-2,cols(d))))|b;
   else;
   a = (1~c)|a;
   b = (1~d)|b;
endif;
retp(a,b,t);
endp;
