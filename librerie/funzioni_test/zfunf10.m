function y=zfunf10(x)
%*****************************************************
% funzione test per la ricerca degli zeri
% y=(1+(1+n)^2)*x-(1-n*x)^2 per n=1,5,10 def. in [0,1]
%*****************************************************
if (nargin==0)
  y(1)=0;
  y(2)=1;
else
  n=10;
  y=(1+(1+n).^2).*x-(1-n.*x).^2;
end
