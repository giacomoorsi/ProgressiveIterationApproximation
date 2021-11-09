function y=zfunf09(x)
%*****************************
% funzione test y=x^3-8 per 
% la ricerca degli zeri [-2,2]
%*****************************
if (nargin==0)
  y(1)=0;
  y(2)=3;
else
  y=x.^3-8;
end
