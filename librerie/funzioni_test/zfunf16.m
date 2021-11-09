function y=zfunf16(x)
%****************************************
% funzione test per la ricerca degli zeri
% y=x^11+4.*x^2-10 def. in [-4,4]
%****************************************
if (nargin==0)
  y(1)=-4;
  y(2)=4;
else
  y=x.^11+4.*x.^2-10;
end
