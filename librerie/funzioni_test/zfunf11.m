function y=zfunf11(x)
%*****************************************
% funzione test per la ricerca degli zeri
% y=x^2-(1-x)^n per n=1,5,10 def. in [0,1]
%*****************************************
if (nargin==0)
  y(1)=0;
  y(2)=1;
else
  n=5;
  y=x.^2-(1-x).^n;
end
