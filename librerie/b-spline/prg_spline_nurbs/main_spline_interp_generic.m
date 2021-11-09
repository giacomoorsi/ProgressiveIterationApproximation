function main_spline_interp_generic(g,n)
% Interpolazione spline di grado generico 
% g   : grado della spline
% n   : n+1 punti di interpolazione
% viene prodotto il grafico della funzione test, della
% spline interpolante e della funzione errore

close all
%setting della funzione test
fun=@runge;
% fun=@ftest1;

y=fun();
a=y(1);
b=y(2);

%punti equispaziati in [a,b]
x=linspace(a,b,n+1);

%campionamento funzione test
y=fun(x);

% disposizione dei nodi in funzione dei punti inseriti
% secondo de Boor per soddisfare le condizioni di Schoenberg-Whitney
t=an_not_a_knot(g,x);

c=an_bspl(g,t,x);
%     disp(c)
%     norm(c,Inf)
%     L=eig(c)
p=c\y';

nxx=201;
xx=linspace(t(g+1),t(end-g),nxx);
% ni=21;
% [nxx,xx]=an_mesh(g,t,ni);

yr=fun(xx);
% yy=an_spl2val(p,t,xx);
yy=bspeval(g,p',t,xx);

figure
hold on;
plot(xx,yy,'b','LineWidth',1.5);
plot(xx,yr,'r','LineWidth',1.5);
plot(x,y,'ko');

figure
hold on;
plot(xx,abs(yr-yy),'g','LineWidth',1.5);

Err=max(abs(yr-yy));
fprintf('Max Absolute Error: %e \n',Err);

