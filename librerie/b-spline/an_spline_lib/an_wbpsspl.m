function [bs,bs1,bs2]=an_wbpsspl(g,t,w,x)
% Calcola i valori delle funzioni base RB-spline (non nulle) e delle loro 
% derivate prima e seconda in un punto x;
% se x e' un vettore torna tre matrici; nella riga i-esima di bs ci sono
% i valori delle funzioni RB-spline valutate in x(i); nella riga i-esima
% di bs1 ci sono i valori delle derivate prime delle B-spline in x(i);
% nella riga i-esima di bs2 ci sono i valori delle derivate seconde delle 
% B-spline in x(i);
% function [bs,bs1,bs2]=an_wbpsspl(g,t,w,x)
% g   --> grado della spline
% t   --> vettore dei nodi partizione estesa
% w   --> vettore dei pesi
% x   --> vettore dei punti di valutazione
% bs  <-- matrice dei valori delle funzioni B-spline nei punti
% bs1 <-- matrice dei valorori delle funzioni derivata prima delle B-spline
%         nei punti
% bs2 <-- matrice dei valorori delle funzioni derivata seconda delle B-spline
%         nei punti

% calcola le funzioni base (non nulle) e le loro derivate prima e seconda in un punto;
% se x e' un vettore torna una matrice; in ogni riga di bs ci sono
% i valori delle funzioni base e in ogni riga di bs1 i valori delle loro
% derivate prima in un punto.
% g --> grado della spline
% t --> vettore dei nodi partizione estesa
% w --> vettore dei pesi
% x --> vettore dei punti di valutazione
% bs  <-- matrice delle funzioni B-spline nei punti
% bs1 <-- matrice delle funzioni derivata prima delle B-spline nei punti
% bs2 <-- matrice delle funzioni derivata seconda delle B-spline nei punti
[bs,bs1,bs2]=an_bpsspl(g,t,x);
mpk=length(t)-g-1;
ll=findspan(mpk-1,g,x,t)+1;
nt=length(t);
m=length(x);
for ii=1:m
 k=ll(ii);
 l=k-g;
%  sbs=0.0;
%  sbs1=0.0;
%  sbs2=0.0;
%  for j=l:k
%      bs(ii,j)=w(j)*bs(ii,j);
%      sbs=sbs+bs(ii,j);
%      bs1(ii,j)=w(j)*bs1(ii,j);
%      sbs1=sbs1+bs1(ii,j);
%      bs2(ii,j)=w(j)*bs2(ii,j);
%      sbs2=sbs2+bs2(ii,j);
%  end
%  for j=l:k
%      bs(ii,j)=bs(ii,j)/sbs;
%      bs1(ii,j)=(bs1(ii,j)-bs(ii,j)*sbs1)/sbs;
%      bs2(ii,j)=(bs2(ii,j)-2*bs1(ii,j)*sbs1-bs(ii,j)*sbs2)/sbs;
%  end
bs(ii,l:k)=w(l:k).*bs(ii,l:k);
sbs=sum(bs(ii,l:k));
bs1(ii,l:k)=w(l:k).*bs1(ii,l:k);
sbs1=sum(bs1(ii,l:k));
bs2(ii,l:k)=w(l:k).*bs2(ii,l:k);
sbs2=sum(bs2(ii,l:k));
bs(ii,l:k)=bs(ii,l:k)/sbs;
bs1(ii,l:k)=(bs1(ii,l:k)-bs(ii,l:k)*sbs1)/sbs;
bs2(ii,l:k)=(bs2(ii,l:k)-2*bs1(ii,l:k)*sbs1-bs(ii,l:k)*sbs2)/sbs;
end
