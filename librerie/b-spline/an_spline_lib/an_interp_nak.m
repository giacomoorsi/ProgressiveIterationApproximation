function [p,t]=an_interp_nak(g,x,y)
% Interpolante spline cubica  (not-a-knot) di (n+1)
% (n>=3) punti (x,y) risolvendo il sistema lineare;
% function [p,t]=an_interp_nak(g,x,y)
% g --> grado della spline di interpolazione
% x --> vettore delle ascisse
% y --> vettore delle ordinate
% p <-- vettore dei coefficienti determinati
% t <-- partizione nodale not-a-knot e nodi
%       aggiuntivi coincidenti
n=length(x)-1;
t(1:4)=x(1)*ones(size(1:4));
t(5:n+1)=x(3:n-1);
t(n+2:n+5)=x(n+1)*ones(size(1:4));

% % *** soluzione non ottimizzata ***
% c=an_bspl(g,t,x);
% % disp(c);
% p=c\y';
% % *********************************

% % *** soluzione ottimizzata ****
p(1)=y(1);
p(n+1)=y(n+1);
c=an_bspl(g,t,x(2:n));
y(2)=y(2)-p(1).*c(1,1);
y(n)=y(n)-p(n+1).*c(n-1,n+1);
c=c(:,2:n);
p(2:n)=c\y(2:n)';
% % *********************************
% figure
% spy(c)
% eig(c)

