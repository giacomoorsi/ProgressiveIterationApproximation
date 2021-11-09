function bs=an_wbspl(g,t,w,x)
% Calcola i valori delle funzioni base RB-spline (non nulle) in un punto x;
% se x e' un vettore torna una matrice; nella riga i-esima di bs ci sono
% i valori delle funzioni B-spline valutate in x(i); 
% function bs=an_wbspl(g,t,w,x)
% g  --> grado della nurbs
% t  --> vettore dei nodi partizione estesa
% w  --> vettore dei pesi
% x  --> vettore dei punti di valutazione
% bs <-- matrice dei valori delle funzioni RB-spline nei punti
bs=an_bspl(g,t,x);
mpk=length(t)-g-1;
ll=findspan(mpk-1,g,x,t)+1;
m=length(x);
for ii=1:m
 k=ll(ii);
 l=k-g;
%  sbs=0.0;
%  for j=l:k
%      bs(ii,j)=w(j)*bs(ii,j);
%      sbs=sbs+bs(ii,j);
%  end
%  for j=l:k
%      bs(ii,j)=bs(ii,j)/sbs;
%  end
 bs(ii,l:k)=w(l:k).*bs(ii,l:k);
 sbs=sum(bs(ii,l:k));
 bs(ii,l:k)=bs(ii,l:k)/sbs;
end

