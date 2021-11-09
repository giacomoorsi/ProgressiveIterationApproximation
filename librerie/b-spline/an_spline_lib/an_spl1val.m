function y=an_spl1val(c,t,x)
% Valutazione della funzione spline definita dai coefficienti
% c e dai nodi t nei punti x mediante combinazione lineare
% delle B-spline; chiama function an_bspl.m per le funzioni base
% function y=an_spl1val(c,t,x)
% c --> coefficienti della spline
% t --> vettore dei nodi partizione estesa
% x --> vettore dei punti di valutazione
% y <-- valore della spline nei punti
c=c(:);
n=length(c);
nt=length(t);
g=nt-n-1;
bs=an_bspl(g,t,x);
y=bs*c;

