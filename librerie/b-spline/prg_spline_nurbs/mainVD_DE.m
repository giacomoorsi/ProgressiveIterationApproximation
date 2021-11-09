function mainVD_DE()
% script per spline su cui sperimentare il degree-elevation 
% g --> grado spline
% mpK --> dimensione spazio spline

close all

g=3;
a=0;
b=1;
% inserimento manuale dei coefficienti
Y=[1 3 -1 3 2 -1 0];
mpK=length(Y);

%nodi aggiuntivi coincidenti
t(1:g)=zeros(1,g);
%nodi aggiuntivi distinti
%t(1:g)=linspace(-0.3,-0.1,g);
if(mpK>2)
 t(g+1:mpK+1)=linspace(0,1,mpK-g+1);
end
%nodi aggiuntivi coincidenti
t(mpK+2:mpK+g+1)=ones(1,g);
%nodi aggiuntivi distinti
%t(n+2:n+g+1)=linspace(1.1,1.3,g);

%ascisse di Greville
gp = aveknt(t,g+1);
figure(1);
hold on;
plot(gp,Y,'k-o');

%punti di valutazione
n=31;
[nm,mesh]=an_mesh(g,t,n);

%valutazione spline
YY=bspeval(g,Y,t,mesh);
plot(mesh,YY,'g','LineWidth',1.5);

flag=1
while (flag==1)
fprintf('partizione estesa:\n');
fprintf('%7.3f ',t);
fprintf('\n');

flag=input('vuoi fare degree elevation? (1-si, 0-no): ');
if flag==1

    %elevamento di 1 grado
    [Ynew,tnew] = bspdegelev(g,Y,t,1);
    g=g+1;
    
    %nuove ascisse di Greville
    gpnew = aveknt(tnew,g+1);
    plot(gpnew,Ynew,'b-o');

    %nuova valutazione spline
    YYnew=bspeval(g,Ynew,tnew,mesh);
    plot(mesh,YYnew,'r','LineWidth',1.5);
    Y=Ynew;
    t=tnew;
end
end

