function y=zfunf08(x)
%*******************************
% funzione test y=x-2^(-x) per 
% la ricerca degli zeri in [0,1]
% sol =
%*******************************
if (nargin==0)
  y(1)=0;
  y(2)=1;
else
  y=x-2.^(-x);
end
