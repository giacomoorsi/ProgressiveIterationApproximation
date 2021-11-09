function y=ftest1(x)
% funzione da approssimare
% definita in [-1,1]
if (nargin==0)
    y(1)=-1.0;
    y(2)=1.0;
else
    y=exp(x) .* sin(2*pi*x);
end