function [p,t]=an_interp_der(g,x,y,y1)
% Interpolante spline cubica derivate agli estremi
% function [p,t]=an_interp_der(g,x,y,y1)
% g  --> grado della spline di interpolazione
% x  --> vettore delle ascisse
% y  --> vettore delle ordinate
% y1 --> vettore dei valori di derivata negli estremi
% p  <-- vettore dei coefficienti determinati
% t  <-- partizione nodale coincidente con i punti x e
%        nodi aggiuntivi coincidenti
n=length(x)-1;
t(1:4)=x(1).*ones(size(1:4));
t(5:n+3)=x(2:n);
t(n+4:n+7)=x(n+1).*ones(size(1:4));

p(1)=y(1);
p(2)=1./3.*y1(1).*(x(2)-x(1))+y(1);
p(n+2)=y(n+1)-1./3.*y1(2).*(x(n+1)-x(n));
p(n+3)=y(n+1);
if (n>1)
  c=an_bspl(g,t,x(2:n));
  y(2)=y(2)-p(2).*c(1,2);
  y(n)=y(n)-p(n+2).*c(n-1,n+2);
  c=c(:,3:n+1);
  p(3:n+1)=c\y(2:n)';
end
