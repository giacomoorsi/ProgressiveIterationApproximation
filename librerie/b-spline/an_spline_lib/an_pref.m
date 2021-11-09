function tref=an_pref(g,t)
%Determina il vettore dei nodi relativo al degree-elevation di una spline
%definita sulla partizione nodale estesa t
%(si assume che t abbia i nodi aggiuntivi coincidenti)
% function tref=an_pref(g,t)
% g    --> grado del polinomio
% t    --> vettore dei nodi partizione estesa
% tref <-- vettore nuovo dei nodi partizione estesa

nt=length(t);

a=t(g+1);
%dalla partizione nodale t costruiamo il vettore nuovo dei nodi
[t_mult,t_single,nt_s]=an_knot_mult(g,t);
tref=zeros(1,nt+nt_s+2);
tref(1:g+2)=t(g+1)*ones(1,g+2);
k=g+2;
for i=1:nt_s
   for j=1:t_mult(i)+1
       k=k+1;
       tref(k)=t_single(i);
   end
end
tref(k+1:k+g+2)=t(nt-g)*ones(1,g+2);

