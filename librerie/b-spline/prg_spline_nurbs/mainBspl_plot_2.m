% viene prodotto il grafico delle funzioni base B-spline di grado g sulla 
% partizione nodale estesa t ...
% e delle loro derivate prime e seconda:
% flag=1 solo B-spline
% flag=2 B-spline e derivate prime
% flag=3 B-spline, derivate prime e seconde
clear all
close all

%definizione spazio spline: grado e nodi partizione estesa
g=3;
%t=[0 0 0 0 1 3 3 3 3];
t = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3];


%flag setting
flag=1;

%numero nodi partizione estesa
nt=length(t);

%dimensione spazio spline
mpk=nt-g-1;

%estremi dell'intervallo di definizione
a=t(g+1);
b=t(nt-g);

%griglia di valutazione uniforme su [a,b]
np=100;
x=linspace(a,b,np);

%Codice per grafico funzioni base B-spline
if flag==1
 y=an_bspl(g,t,x);
 figure(1);
 hold on;
 for i=1:mpk
   plot(x,y(:,i),fcol(i),'LineWidth',1.5);
 end
end

% %Codice per grafico B-spline e derivate prima
if flag==2
 [y,y1]=an_bpspl(g,t,x);
 figure(1);
 hold on;
 for i=1:mpk
   plot(x,y(:,i),fcol(i),'LineWidth',1.5);
 end
 figure(2);
 hold on;
 for i=1:mpk
   plot(x,y1(:,i),fcol(i),'LineWidth',1.5);
 end
end

% Codice per grafico B-spline, derivate prima e seconda
if flag==3
 [y,y1,y2]=an_bpsspl(g,t,x);
 figure(1);
 hold on;
 for i=1:mpk
   plot(x,y(:,i),fcol(i),'LineWidth',1.5);
 end
 figure(2);
 hold on;
 for i=1:mpk
   plot(x,y1(:,i),fcol(i),'LineWidth',1.5);
 end
 figure(3);
 hold on;
 for i=1:mpk
   plot(x,y2(:,i),fcol(i),'LineWidth',1.5);
 end
end

%disegno partizione nodale
plot(t,zeros(1,length(t)),'k*-','LineWidth',1.5);

function ic=fcol(k)
%definizione colori
col=['k','r','g','b','m','c','y'];
ic=col(mod(k-1,7)+1);
end