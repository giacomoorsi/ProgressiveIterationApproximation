function y=zfunf05(x)
%**********************************
% funzione test y=x^3-3*x+2 per 
% la ricerca degli zeri in [-2.5,2]
% sol1 = -2, sol2 = 1
%**********************************
if (nargin==0)
  y(1)=-2.5;
  y(2)=2;
else
  y=x.^3-3.*x+2;
end
