function bs=an_bspl(g,t,x)
% Calcola i valori delle funzioni base B-spline (non nulle) in un punto x;
% se x e' un vettore torna una matrice; nella riga i-esima di bs ci sono
% i valori delle funzioni B-spline valutate in x(i); 
% function bs=an_bspl(g,t,x)
% g  --> grado della spline
% t  --> vettore dei nodi partizione estesa
% x  --> vettore dei punti di valutazione
% bs <-- matrice dei valori delle funzioni B-spline nei punti
mpk=length(t)-g-1;
iv = findspan(mpk-1,g,x,t);
bb = basisfun(iv,x,g,t);
gp1=1:g+1;
for j=1:length(x)
   jj=iv(j)+2-flip(gp1);
   bs(j,jj)=bb(j,gp1);
%    bs(j,iv(j)-g+1:iv(j)+1)=bb(j,1:g+1);
end
