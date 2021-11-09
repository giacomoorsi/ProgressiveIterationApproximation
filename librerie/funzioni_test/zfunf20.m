function y=zfunf20(x)
%****************************************
% funzione test per la ricerca degli zeri
% y=cosh(x)+cos(x)-gamma def. in [0,1]
% con gamma=1,2,3
%****************************************
if (nargin==0)
  y(1)=-2;
  y(2)=2;
else
  gamma=2;
  y=cosh(x)+cos(x)-gamma;
end
