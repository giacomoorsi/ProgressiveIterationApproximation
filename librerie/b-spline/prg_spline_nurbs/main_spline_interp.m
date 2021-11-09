function main_spline_interp(n,flag)
% Interpolazione spline cubica 
% n   : n+1 punti di interpolazione
% flag: scelta del metodo/tipo di interpolazione
% 1: cubica not-a-knot
% 2: cubica con derivate agli estremi
% 3: cubica periodica
% 4: cubica naturale
% 5: cubica di Hermite
% viene prodotto il grafico della funzione test, della spline interpolante
% e della funzione errore

close all

%setting della funzione test
fun=@runge;
d1fun=@d1runge;

%acquisisce gli estremi dell'intervallo di definizione
y=fun();
a=y(1);
b=y(2);

%punti equispaziati in [a,b]
x=linspace(a,b,n+1);

%campionamento funzione test
y=fun(x);

%grado della funzione spline
g=3;

%metodo di interpolazione
switch flag
    case 1
%condizioni 'Not a Knot'
       [p,t]=an_interp_nak(g,x,y);
    case 2
%condizioni 'derivate agli estremi'
       y1(1)=d1fun(x(1));
       y1(2)=d1fun(x(n+1));
       [p,t]=an_interp_der(g,x,y,y1);
    case 3
%condizioni 'periodiche'
       [p,t]=an_interp_per(g,x,y);
    case 4
%condizioni 'naturali'
       [p,t]=an_interp_nat(g,x,y);
    case 5
%derivate campionate dalla funzione derivata
       dy=d1fun(x);
%cubica di Hermite C^1
       [p,t]=an_interp_SH3(x,y,dy);
end

ni=21;
[nxx,xx]=an_mesh(g,t,ni);
yr=fun(xx);
% yy=an_spl2val(p,t,xx);
yy=bspeval(g,p,t,xx);

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

