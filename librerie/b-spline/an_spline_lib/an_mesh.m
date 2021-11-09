function [nm,mesh]=an_mesh(g,t,n)
%Genera una mesh opportuna dell'intervallo di definizione rispettando
%i nodi e le loro molteplicita'
%function [nm,mesh]=an_mesh(g,t,n)
%g --> grado della spline
%t --> partizione nodale estesa
%n --> numero di punti della mesh per intervallo nodale
%nm    <-- numero punti mesh generati
%mesh  <-- mesh generata
if (n<2)
  n=2;
end
nt=length(t);
j=1;
for i=g+1:nt-g-1
  if (t(i+1)>t(i))
   nm=j+n-1;
   mesh(j:nm)=linspace(t(i),t(i+1),n);
%modifica affinche' i nodi della partzione estesa 
%siano duplicati nella mesh
  j=nm;
%    j=nm+1;
  end
end