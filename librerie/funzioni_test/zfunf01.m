function y=zfunf01(x)
%****************************************
% funzione test per la ricerca degli zeri
% y=(4*x-7)/(x-2) def. in [1,1.9]
% sol = 7/4 = 1.75
%****************************************
if (nargin==0)
    y(1)=1.0;
    y(2)=1.9;
else
    y=(4.*x-7)./(x-2);
end
