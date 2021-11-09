function y=zfunf19(x)
%****************************************
% funzione test per la ricerca degli zeri
% y=(x-eps) def. in [0,1]
% con eps=2.220446049250313e-16 predefinito in Matlab
%****************************************
if (nargin==0)
  y(1)=0;
  y(2)=1;
else
  y=(x-eps);
end
