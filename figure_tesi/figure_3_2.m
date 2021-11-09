
clear all
close all

fun = 'epitrochoid';
a = 0;
b = 6*pi;


g = 3; % grado della spline
grafici = 3; % 3 per mostrare tutti i grafici
n_q = 100; % numero punti di approssimazione
n = 24; % numero punti di controllo
equidistanti = 1; % 1 per estrarre punti equidistanti dai punti campionati
knot_partition = 3;

q_tt = linspace(a, b, n_q);
[q_x, q_y] = feval(fun, q_tt);



% prendo n+1 punti tra quelli campionati (saranno i punti iniziali LSPIA)
indici = floor(linspace(1, length(q_tt), n+1));

x = q_x(indici);
y = q_y(indici);

tt = q_tt(indici);

tol1 = 1e-7; % tolleranza nei punti
tol2 = 1e-15; % avanzamento minimo

adjust_axis = 1; % 1 per gli svg

a = min(tt);
b = max(tt);


chord_parametrization = 1;

% LSPIA per B-Spline
%
% parametri da impostare prima della chiamata a lspia_bspline_body:
% q_x, q_y              <- vettori punti di approssimazione
% q_tt                  <- parametri dei punti di approssimazione
% x, y                  <- vettore punti 
% tt                    <- parametri dei punti x, y
% tol1                  <- tolleranza sull'errore di approssimazione
% tol2                  <- tolleranza sull'avanzamento (diff. errori al
%                       passo k e al passo k-1)
% grafici               <- 1 per conclusivo, 3 per grafici su più passi
% knot_partition:
%                       1. per prendere i parametri dei punti
%                       2. per punti uniformi in tt(1)...tt(end)
%                       3. per punti  Whitney–Schoenberg
% g                     <- grado della B-Spline
% adjust_axis           <- 1 per svg



if (~exist('adjust_axis'))
    adjust_axis = 0;
end


p = [x;y];
a = min(tt);
b = max(tt);

% iterazioni da mostrare nel subplot (in ultima posizione è sempre mostrata
% l'ultima iterazione)
step_to_show = [0, 1, 3, 5, 10, 30, 50, 80, 100, 200];

k = length(p)-g-1;
if (knot_partition == 1)
    t = [repmat(a,1,g+1), tt(3:end-2), repmat(b,1,g+1)];
elseif (knot_partition == 2)
    t = [repmat(a,1,g), linspace(a, b, k+2), repmat(b,1,g)];
elseif (knot_partition == 3)
    t = an_not_a_knot(g,tt);
end
bs = an_bspl(g, t, q_tt);


btb = bs' * bs;

t_val = linspace(tt(1), tt(end), 3000);
bs_val_err = an_bspl(g, t, t_val);

[x_real_err, y_real_err] = feval(fun, t_val);
x_real_err = x_real_err';
y_real_err = y_real_err';


    
tabella = [] %


% Peso esatto
% autoval = eig(btb);
% mu = 2/(min(autoval)+max(autoval));

% Peso ottimizzato
C = max(sum(btb, 2)); % prendo il massimo della somma di ogni riga
mu = 2/C; % definisco mu come da approssimazione

err1 = tol1 + 1;
avanzamento = tol2 + 1;
k = 0;

if(grafici > 0)
    n_val = 5000;
    p_val = linspace(a, b, n_val);
    bs_val = an_bspl(g, t, p_val);
end
 if (grafici==3 || any(ismember(step_to_show, k)))
        n_grafico = 1;
        %subplot(3,2,n_grafico)
        figure(1)
        plot((bs_val*p(1,:)')', (bs_val*p(2,:)')','r');
        hold on
        plot(q_x, q_y, 'b.');
        plot(p(1,:), p(2,:),'g+');
        axis([-2.8 2.8 -2.8 2.8])
      %  title(strcat('k: ', num2str(k)))
        
        if (adjust_axis)
            axis ij
            axis equal
        end
    end

% itero finché l'errore di approssimazione è superiore alla tolleranza e
% finché il passo di avanzamento è anch'esso superiore alla tolleranza

while(err1 > tol1 && avanzamento > tol2)
    k = k+1;
    delta = [q_x - (bs*p(1,:)')'; q_y - (bs*p(2,:)')'];
    for(i=1:(n+1))
        adj(1, i) = mu * sum(bs(:, i)'.*delta(1, :));
        adj(2, i) = mu * sum(bs(:, i)'.*delta(2, :));
    end
    
    p = p + adj;
    
    for j=1:length(delta)
       norma(j) = norm(delta(:,j)) ;
    end
    tmp = max(norma);
    avanzamento = abs(err1 - tmp);
    err1 = tmp;
    
    if(k<10 || mod(k-1, 5)==0)
        fprintf("Errore al passo %1d: %5.5e\n", (k-1), err1);
        fprintf("Avanzamento al passo %1d: %5.5e\n", (k-1), avanzamento);
    end
    
    if (grafici==3 && any(ismember(step_to_show, k)))
        n_grafico = n_grafico + 1;
        %subplot(3,2,n_grafico)
        figure(1)
        hold off
        plot((bs_val*p(1,:)')', (bs_val*p(2,:)')','r');
        hold on
        plot(q_x, q_y, 'b.');
        plot(p(1,:), p(2,:),'g+');
        axis([-2.8 2.8 -2.8 2.8])
       % title(strcat('k: ', num2str(k)))
        
    end
    
    x_interp_err = bs_val_err*p(1,:)';
    y_interp_err = bs_val_err*p(2,:)';
    tabella(k+1, 1) = k-1;
    tabella(k+1, 2) = err1; 
    tabella(k+1, 3) = max_distance([x_interp_err y_interp_err], [x_real_err y_real_err]);
    
end
if(avanzamento <= tol2)
   fprintf("Interrotto per tolleranza raggiunta su avanzamento (%3.0e)\n", tol2); 
end


if (grafici > 0)
    figure(1)
    hold off
    plot(q_x, q_y, 'b.', 'DisplayName', 'Punti di approssimazione');
    hold on
    plot((bs_val*p(1,:)')', (bs_val*p(2,:)')','r', 'DisplayName', 'Curva approssimata');
    plot(p(1,:), p(2,:), 'g+', 'DisplayName', 'Punti di controllo');
    %legend
    axis([-2.8 2.8 -2.8 2.8])
    if (adjust_axis)
        axis ij
        axis equal
    end
    
    %plottitle = strcat('n-Q: ', num2str(n_q), ' - n: ', num2str(n), ' - k: ', num2str(k), ' - tol1: ', num2str(tol1), ' - tol2: ', num2str(tol2));
    %title(plottitle)
    
    fprintf("Passo %1d. Errore: %5.5e.\n", k, err1);
end


function [x,y] = epitrochoid(t)
% definita tra 0 e 6*pi
a=2;
b=3;
c=1;
x = (a-b)*cos(t)-c*cos((a/b+1)*t);
y = (a-b)*sin(t)-c*sin((a/b+1)*t);
end


