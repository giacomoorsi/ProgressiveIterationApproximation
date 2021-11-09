function y=zfunf15(x)
%*******************************************************
% funzione test per la ricerca degli zeri
% y=tanh(x-1.0) def. in [-1,3]
%*******************************************************
if (nargin==0)
  y(1)=-1;
  y(2)=3;
else
  y=tanh(x-1.0);
end
