function y=zfunf04(x)
%*****************************
% funzione test y=x^n-2 per la 
% ricerca degli zeri in [0,2]
% n=2: sol = 2^(1/2)
% n=3: sol = 2^(1/3)
% n=4: sol = 2^(1/4)
%*****************************
if (nargin==0)
  y(1)=0;
  y(2)=2;
else
  n=2;
  y=x.^n-2;
end
