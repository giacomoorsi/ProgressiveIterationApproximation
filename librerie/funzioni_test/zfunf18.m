function y=zfunf18(x)
%****************************************
% funzione test per la ricerca degli zeri
% y=(20*x-135)/(x-7) def. in [6,7]
% sol1=6.75
%****************************************
if (nargin==0)
    y(1)=6;
    y(2)=7;
else
%     y=(20.*x-135)./(x-7);
    y=20.*x-135;
end
