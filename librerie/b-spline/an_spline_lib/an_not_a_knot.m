function t=an_not_a_knot(g,x)
%Dispone i nodi rispetto ai punti x con il
%metodo proposto da de Boor; nel caso g=3
%questo corrisponde al 'not-a-knot'
%function t=an_not_a_knot(g,x)
%g --> grado della spline
%x --> punti di interpolazione
%t <-- partizione estesa dei nodi secondo de Boor
m=g+1;
mpk=length(x);
t=zeros(1,mpk+m);
t(1:m)=x(1)*ones(1,m);
for i=m+1:mpk
 t(i)=sum(x(i+1-m:i-1))/g;
end
t(mpk+1:mpk+m)=x(mpk)*ones(1,m);

