function y=runge(x)
% funzione di runge
% y=1./(1.+x.^2)
% x in [-5,5]
if (nargin==0)
 y(1)=-5.0;
 y(2)=5.0;
else
 y=1./(1.+x.^2);
end
