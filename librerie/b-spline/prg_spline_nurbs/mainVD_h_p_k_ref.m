function mainVD_h_p_k_ref()
% Valutazione e grafico di spline VD 

clear all
close all

%Colori
col=['k','r','g','b','m','c','y','k','r','g','b','m','c','y'];

a=0;
b=1;
% inserimento manuale dei coefficienti
hold on;
g=1;
t=[0 0 1 1];
Y=[1 1];
mpK=length(Y);

figure(1);
hold on;
n=31;
[nm,mesh]=an_mesh(g,t,n);

iv = findspan(mpK,g,mesh,t);
iv(end)=iv(end-1);
yy = basisfun(iv,mesh,g,t);
for j=1:nm
   y(j,iv(j)-g+1:iv(j)+1)=yy(j,1:g+1);
end

figure(1);
hold on;
for i=1:mpK
    plot(mesh,y(:,i),col(i),'LineWidth',1.5);
end

nstep=1;
for kk=1:nstep

%h-refinement
    [tref,tc]=an_href(g,t);
    [Y,t]=bspkntins(g,Y,t,tc);
   
%p-refinement
    [Y,t] = bspdegelev(g,Y,t,1);
    g=g+1;
    
% %h-refinement
%     [tref,tc]=an_href(g,t);
%     [Y,t]=bspkntins(g,Y,t,tc);

end

n=31;
[nm,mesh]=an_mesh(g,t,n);

mpK=length(Y);
 
y  = an_bspl(g,t,mesh);

figure(2);
hold on;
for i=1:mpK
 plot(mesh,y(:,i),col(i),'LineWidth',1.5);
end

