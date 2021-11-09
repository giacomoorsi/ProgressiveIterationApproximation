function y=zfunf02(x)
%****************************************
% funzione test per la ricerca degli zeri
% y=1+3/x^2-4/x^3 def. in [0.5,6]
% sol = 1.0
%****************************************
if (nargin==0)
    y(1)=0.5;
    y(2)=6;
else
    y=1+3./x.^2-4./x.^3;
end
