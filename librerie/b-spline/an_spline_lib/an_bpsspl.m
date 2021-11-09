function [bs,bs1,bs2]=an_bpsspl(g,t,x)
% Calcola i valori delle funzioni base B-spline (non nulle) e delle loro 
% derivate prima e seconda in un punto x;
% se x e' un vettore torna tre matrici; nella riga i-esima di bs ci sono
% i valori delle funzioni B-spline valutate in x(i); nella riga i-esima
% di bs1 ci sono i valori delle derivate prime delle B-spline in x(i);
% nella riga i-esima di bs2 ci sono i valori delle derivate seconde delle 
% B-spline in x(i);
% function [bs,bs1,bs2]=an_bpsspl(g,t,x)
% g   --> grado della spline
% t   --> vettore dei nodi partizione estesa
% x   --> vettore dei punti di valutazione
% bs  <-- matrice dei valori delle funzioni B-spline nei punti
% bs1 <-- matrice dei valorori delle funzioni derivata prima delle B-spline
%         nei punti
% bs2 <-- matrice dei valorori delle funzioni derivata seconda delle B-spline
%         nei punti
mpk=length(t)-g-1;
iv = findspan(mpk-1,g,x,t);
bb = basisfunder(iv,g,x,t,2);
gp1=1:g+1;
for j=1:length(x)
   jj=iv(j)+2-flip(gp1);
   bs(j,jj)=bb(j,1,gp1);
   bs1(j,jj)=bb(j,2,gp1);
   bs2(j,jj)=bb(j,3,gp1);
end