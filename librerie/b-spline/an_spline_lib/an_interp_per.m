function [p,t]=an_interp_per(g,x,y)
% Interpolante spline cubica periodica di (n+1)
% (n>=3) punti (x,y) risolvendo il sistema lineare;
% function [p,t]=an_interp_per(g,x,y)
% g --> grado della spline di interpolazione
% x --> vettore delle ascisse
% y --> vettore delle ordinate
% p <-- vettore dei coefficienti determinati
% t <-- partizione nodale coincidente con i punti x e nodi
%       aggiuntivi distinti per periodica
n=length(x)-1;
%punti equidistanti in [a,b]
t(4:n+4)=x;
t(3)=t(4)-(t(n+4)-t(n+3));
t(2)=t(4)-(t(n+4)-t(n+2));
t(1)=t(4)-(t(n+4)-t(n+1));
t(n+5)=t(n+4)+(t(5)-t(4));
t(n+6)=t(n+4)+(t(6)-t(4));
t(n+7)=t(n+4)+(t(7)-t(4));

c=an_bspl(g,t,x(1:n));
c(1,n+1)=c(1,n+1)+c(1,1);
c(n,2)=c(n,2)+c(n,n+2);
c=c(:,2:n+1);
p(2:n+1)=c\y(1:n)';
p(1)=p(n+1);
p(n+2)=p(2);
p(n+3)=p(3);