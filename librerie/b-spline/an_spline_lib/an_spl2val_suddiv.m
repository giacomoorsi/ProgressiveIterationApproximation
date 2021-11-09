function [c1,c2]=an_spl2val_suddiv(c,t,x)
% Suddivisione della funzioni spline definita dai coefficienti
% c e dai nodi t nel punto x interno all'intervallo, mediante l'algoritmo 
% di de Boor sui coefficienti
% function y=an_spl2val_suddiv(c,t,x)
% c  --> coefficienti della spline
% t  --> vettore dei nodi partizione estesa
% x  --> vettore dei punti di valutazione
% c1 <-- coefficienti spline primo intervallo
% c2 <-- coefficienti spline secondo intervallo
n=length(c);
nt=length(t);
g=nt-n-1;
m=length(x);
l=gc_findint(g,t,x);
for k=1:m
 w=c;
 y1(1)=w(l(k)-g);
 y2(g+1)=w(l(k));
 for j=1:g
  for i=l(k):-1:l(k)-g+j
   d1=x(k)-t(i);
   d2=t(i+g-j+1)-x(k);
   w(i)=(d1.*w(i)+d2.*w(i-1))./(d1+d2);
   c1(i)=w(i);
  end
  c2(g-j+1)=w(l(k));
 end
end

