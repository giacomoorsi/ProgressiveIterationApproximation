function y=zfunf13(x)
%*******************************************************
% funzione test per la ricerca degli zeri
% y=1+2*(x*exp(-n)-exp(-n*x)) per n=1,5,10 def. in [0,1]
%*******************************************************
if (nargin==0)
  y(1)=0;
  y(2)=1;
else
  n=5;
  y=1+2.*(x.*exp(-n)-exp(-n.*x));
end
