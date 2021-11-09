function [c,t]=an_interp_SH3(qx,qy,dy)
% Interpolante spline cubica  C^1 di Hermite
% function [c,t]=an_interp_SH3(qx,qy,dy)
% qx --> vettore delle ascisse
% qy --> vettore delle ordinate
% dy --> vettore dei valori derivati
% c  <-- vettore dei coefficienti determinati
% t  <-- partizione nodale e nodi con aggiuntivi coincidenti
n=length(qx);
g=3;
nt=2*n+g+1;
t(1:g+1)=qx(1)*ones(size(1,g+1));
k=2;
for i=g+2:2:nt-g-1
    t(i)=qx(k);
    t(i+1)=qx(k);
    k=k+1;
end
t(nt-g:nt)=qx(n)*ones(size(1,g+1));

c(1)=qy(1);
c(2)=qy(1)+1/3*(qx(2)-qx(1))*dy(1);
for i=2:n-1
    c(2*i-1)=qy(i)-1/3*(qx(i)-qx(i-1))*dy(i);
    c(2*i)=qy(i)+1/3*(qx(i+1)-qx(i))*dy(i);
end
c(2*n-1)=qy(n)-1/3*(qx(n)-qx(n-1))*dy(n);
c(2*n)=qy(n);

end
