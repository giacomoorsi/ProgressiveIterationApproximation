function [tref,tins]=an_href(g,t)
%Genera il vettore dei nodi (partizione estesa) raffinamento del vettore 
%dei nodi t 
% function [tref,tins]=an_href(g,t)
% g    --> grado della spline
% t    --> vettore dei nodi partizione estesa
% tref <-- vettore dei nodi raffinato partizione estesa
% tins <-- vettore di nodi inseriti
nt=length(t);

a=t(g+1);
%dalla partizione nodale t costruiamo il vettore dei nodi raffinati
[t_mult,t_single,nt_s]=an_knot_mult(g,t);
tref=zeros(1,nt+nt_s+1);
tref(1:g+1)=t(1:g+1);
h=0;
k=g+1;
for i=1:nt_s
   h=h+1;
   tins(h)=0.5*(a+t_single(i));
   k=k+1;
   tref(k)=tins(h);
   a=t_single(i);
   for j=1:t_mult(i)
       k=k+1;
       tref(k)=t_single(i);
   end
end
tins(h+1)=0.5*(a+t(nt-g));
tref(k+1)=tins(h+1);
tref(k+2:k+g+2)=t(nt-g:nt);
end
