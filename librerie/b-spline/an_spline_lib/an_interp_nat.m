function [p,t]=an_interp_nat(g,x,y)
% Interpolante spline cubica naturale di (n+1)
% (n>=3) punti (x,y) risolvendo il sistema lineare;
% function [p,t]=an_interp_nat(g,x,y)
% g --> grado della spline di interpolazione
% x --> vettore delle ascisse
% y --> vettore delle ordinate
% p <-- vettore dei coefficienti determinati
% t <-- partizione nodale coincidente con i punti x e
%       nodi aggiuntivi coincidenti
n=length(x)-1;
%punti equidistanti in [a,b]
t(1:4)=x(1).*ones(size(1:4));
t(5:n+3)=x(2:n);
t(n+4:n+7)=x(n+1).*ones(size(1:4));
p(1)=y(1);
p(n+3)=y(n+1);
c=zeros(n+1);
if (n>1)
  cc=an_bspl(g,t,x(2:n));
  c(2:n,:)=cc(:,2:n+2); 
  c(1,1)=2*x(1)-x(2)-x(3);
  c(1,2)=x(2)-x(1);
  c(1,3:n+1)=0;
  c(n+1,n)=x(n+1)-x(n);
  c(n+1,n+1)=x(n-1)+x(n)-2*x(n+1);
  y(1)=p(1)*(x(1)-x(3));
  y(n+1)=p(n+3)*(x(n-1)-x(n+1));
%  disp(c);
  q=c\y';
  p(2:n+2)=q(1:n+1);
end
