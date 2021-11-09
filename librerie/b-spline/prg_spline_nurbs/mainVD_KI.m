function mainVD_KI()
% Valutazione e grafico di spline, rappresentazione 
% grafica dei coefficienti e KI

close all

%grado spline
g=3;
%intervallo di definizione
a=0;
b=1;
% inserimento manuale dei coefficienti
Y=[1 3 -1 3 2 -1 0];
%dimensione spazio spline
mpK=length(Y);

%nodi aggiuntivi coincidenti
t(1:g)=zeros(1,g);
% %nodi aggiuntivi distinti
% t(1:g)=linspace(-0.3,-0.1,g);
if(mpK>2)
 t(g+1:mpK+1)=linspace(0,1,mpK-g+1);
end
%nodi aggiuntivi coincidenti
t(mpK+2:mpK+g+1)=ones(1,g);
% %nodi aggiuntivi distinti
% t(mpK+2:mpK+g+1)=linspace(1.1,1.3,g);

%ascisse di Greville
gp = aveknt(t,g+1);

figure(1);
hold on
plot(gp,Y,'k-o');

%punti di valutazione
n=31;
[nm,mesh]=an_mesh(g,t,n);

%valutazione funzione spline
YY=bspeval(g,Y,t,mesh);
plot(mesh,YY,'g','LineWidth',1.5);

flag=1;
while (flag==1)

fprintf('partizione estesa:\n');
fprintf('%7.3f ',t);
fprintf('\n');

flag=input('vuoi inserire? (1-si, 0-no): ');
if flag==1
    tc=input('nodo/i da inserire: ');

    [Ynew,tnew]=bspkntins(g,Y,t,tc);
    fprintf('partizione estesa \n');
    fprintf('%7.3f ',tnew);
    fprintf('\n');

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
