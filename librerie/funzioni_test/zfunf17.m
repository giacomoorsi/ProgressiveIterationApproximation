function y=zfunf17(x)
%****************************************
% funzione test per la ricerca degli zeri
% y=sin(x)-x def. in [-1,1]
%****************************************
if (nargin==0)
  y(1)=-1;
  y(2)=1;
else
  y=sin(x)-x;
end
