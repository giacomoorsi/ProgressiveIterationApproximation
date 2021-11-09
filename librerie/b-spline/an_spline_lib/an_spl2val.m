function y=an_spl2val(c,t,x)
% Valutazione della funzione spline definita dai coefficienti
% c e dai nodi t nei punti x mediante l'algoritmo di de Boor
% sui coefficienti
% function y=an_spl2val(c,t,x)
% c --> coefficienti della spline
% t --> vettore dei nodi partizione estesa
% x --> vettore dei punti di valutazione
% y <-- valore della spline nei punti
n=length(c);
nt=length(t);
g=nt-n-1;
m=length(x);
l=findspan(n-1,g,x,t)+1;
% l=gc_findint(g,t,x);
for k=1:m
 w=c;
 for j=1:g
  for i=l(k):-1:l(k)-g+j
   d1=x(k)-t(i);
   d2=t(i+g-j+1)-x(k);
   w(i)=(d1.*w(i)+d2.*w(i-1))./(d1+d2);
  end
 end
 y(k)=w(l(k));
end

