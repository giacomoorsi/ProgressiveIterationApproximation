function y=zfunf06(x)
%******************************
% funzione test y=1/x-2 per la 
% ricerca dedgli zeri in [0,4]
% sol =
%******************************
if (nargin==0)
  y(1)=0.1;
  y(2)=4;
else
  y=1./x-2;
end
