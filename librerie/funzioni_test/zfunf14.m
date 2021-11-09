function y=zfunf14(x)
%*******************************************************
% funzione test per la ricerca degli zeri
% y=exp(x)-x^2 def. in [-2,2]
% sol = ?0.7034674
%*******************************************************
if (nargin==0)
  y(1)=-2;
  y(2)=2;
else
  y=exp(x)-x.^2;
end
