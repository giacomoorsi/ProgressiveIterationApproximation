function main_spline_lsq(n,N)
% Approssimante spline cubica (spline di dim. n) nel senso dei
% minimi quadrati discreti di una funzione test;
% viene prodotto il grafico della funzione test, della
% spline approsimante e della funzione errore
% vengono campionati N punti con (N>n>=4) della funzione test;
% viene risolto il sistema delle equazioni normali;
% viene prodotto il grafico della spline approssimante;
% input:
% n --> dimensione spazio spline cubico
% N --> numero di punti di approssimazione 

close all

%setting della funzione test
% fun=@frunge;
fun=@ftest1;
y=fun();
a=y(1);
b=y(2);

g=3;
%punti equidistanti in [a,b]
x=linspace(a,b,N);
y=fun(x);

%nodi aggiuntivi iniziali coincidenti con x(1)
t(1:g)=x(1).*ones(size(1:g));

% disposizione automatica dei knot interni rispetto ai punti di
% approssimazione prendendoli equispaziati come indici
% sul set di punti x dato 
% step=(N-1)/(n-g);
% for i=0:n-g
%  t(i+g+1)=x(round(i*step)+1);
% end

% disposizione automatica dei knot interni rispetto ai punti di
% approssimazione prendendoli equispaziati
t(g+1:n+1)=linspace(a,b,n-g+1);

%nodi aggiuntivi finali coincidenti con x(N)
t(n+2:n+g+1)=x(N).*ones(size(1:g));

c=an_bspl(g,t,x);
p=c\y';

figure(1);
hold on;
yy=an_spl2val(p,t,x);
res=sum((yy-y).^2);
fprintf('Residuo LSQ discreto: %e\n',res);

ni=51;
[nxx,xx]=an_mesh(g,t,ni);
yf=fun(xx);
plot(xx,yf,'g-','LineWidth',1.5);
yy=an_spl2val(p,t,xx);
plot(x,y,'+r');
plot(xx,yy,'r-','LineWidth',1.5);

figure(2);
hold on;
plot(xx,abs(yf-yy),'g','LineWidth',1.5);

w=yf-yy;
Inf_err=max(abs(w));
fprintf('Errore in Norma Inf.: %e\n',Inf_err); 
