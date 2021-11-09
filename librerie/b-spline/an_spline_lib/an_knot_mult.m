function [t_mult,t_single,nt_s]=an_knot_mult(g,t)
%Analizza la partizione nodale e restituisce informazioni
%sui nodi multipli e singoli
% function [t_mult,t_single,nt_s]=an_knot_mult(g,t)
% g --> grado spline
% t --> partizione estesa dei nodi
% t_mult   <-- molteplicita' dei nodi
% t_single <-- break point
% nt_s     <-- numero dei break point
nt=length(t);
t_mult=[];
t_single=[];
nt_s=0;
m=1;
for i=g+2:nt-g-1
   if (t(i) < t(i+1))
      nt_s=nt_s+1;
      t_mult(nt_s)=m;
      t_single(nt_s)=t(i);
      m=1;
   else
       m=m+1;
   end
end
