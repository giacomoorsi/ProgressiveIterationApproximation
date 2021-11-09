function [bs,bs1]=an_wbpspl(g,t,w,x)
% Calcola i valori delle funzioni base RB-spline (non nulle) e delle loro 
% derivate prima in un punto x;
% se x e' un vettore torna due matrici; nella riga i-esima di bs ci sono
% i valori delle funzioni B-spline valutate in x(i); nella riga i-esima
% di bs1 ci sono i valori delle derivate prime delle B-spline in x(i);
% function [bs,bs1]=an_wbpspl(g,t,w,x)
% g   --> grado della spline
% t   --> vettore dei nodi partizione estesa
% w   --> vettore dei pesi
% x   --> vettore dei punti di valutazione
% bs  <-- matrice dei valori delle funzioni RB-spline nei punti
% bs1 <-- matrice dei valorori delle funzioni derivata prima delle B-spline
%         nei punti
[bs,bs1]=an_bpspl(g,t,x);
mpk=length(t)-g-1;
ll=findspan(mpk-1,g,x,t)+1;
m=length(x);
for ii=1:m
 k=ll(ii);
 l=k-g;
%  sbs=0.0;
%  sbs1=0.0;
%  for j=l:k
%      bs(ii,j)=w(j)*bs(ii,j);
%      sbs=sbs+bs(ii,j);
%      bs1(ii,j)=w(j)*bs1(ii,j);
%      sbs1=sbs1+bs1(ii,j);
%  end
%  for j=l:k
%      bs(ii,j)=bs(ii,j)/sbs;
%      bs1(ii,j)=(bs1(ii,j)-bs(ii,j)*sbs1)/sbs;
%  end
bs(ii,l:k)=w(l:k).*bs(ii,l:k);
sbs=sum(bs(ii,l:k));
bs1(ii,l:k)=w(l:k).*bs1(ii,l:k);
sbs1=sum(bs1(ii,l:k));
bs(ii,l:k)=bs(ii,l:k)/sbs;
bs1(ii,l:k)=(bs1(ii,l:k)-bs(ii,l:k)*sbs1)/sbs;
end

