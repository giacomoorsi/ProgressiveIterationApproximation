function bs=an_bernst(g,x)
% Calcola i valori delle funzioni base di Bernstein definite in [0,1] 
% in un punto x;
% se x e' un vettore torna una matrice; nella riga i-esima di bs ci sono
% i valori delle funzioni base di Bernstein valutate in x(i);
% function bs=an_bernst(g,x)
% g  --> grado del polinomio
% x  --> vettore di punti
% bs <-- matrice dei valori delle funzioni base di Bernstein nei punti

m=length(x);
n=g+1;
bs=zeros(m,n);
for k=1:m
 l=n;
 bs(k,l)=1.0;
 d1=x(k);
 d2=1.0-x(k);
 for i=1:g
   temp=0.0;
   for j=l:n
     bs(k,j-1)=d2.*bs(k,j)+temp;
     temp=d1.*bs(k,j);
   end
   bs(k,n)=temp;
   l=l-1;
  end
end

