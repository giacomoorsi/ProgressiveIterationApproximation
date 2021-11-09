
close all
clear all
clc


n = 18;
% funzione
fun = 'helix';
use_optimal_weight = 1;
tol = 1e-16;
grafici = 1;

for k=0:n
    punti(k+1) = k * (pi/3);
    [x(k+1), y(k+1), z(k+1)] = feval(fun, punti(k+1));
end
a = min(punti);
b = max(punti);

k = 0;
g = 3;

p = [x;y;z];

k = length(p)-g-1;
t = an_not_a_knot(g,punti);

bs = an_bspl(g,t,punti);

if (use_optimal_weight == 1)
    autoval = eig(bs);
    min_autoval = min(abs(autoval));
    w = 2 / (1 + min_autoval);
else
    w = 1;
end

d = [x - (bs*p(1,:)')'; y-(bs*p(2,:)')'; z-(bs*p(3,:)')'];

p_val = linspace(a, b, 100);
bs_val = an_bspl(g, t, p_val);


for j=1:length(d)
    norma(j) = norm(d(:,j));
end
err = max(norma);
k = 0;
fprintf("Errore al passo %1d: %5.5e\n", k, err);

if(grafici > 0)
    [x_real, y_real, z_real] = feval(fun, p_val);
    x_interp = bs_val*p(1,:)';
    y_interp = bs_val*p(2,:)';
    z_interp = bs_val*p(3,:)';
    
    close all
    
    %figure(1)
    %plot3(x_real, y_real, z_real, 'r-', 'DisplayName', 'Funzione interpolata');
    %hold on
    %plot3(x_interp, y_interp,  z_interp,'b-', 'DisplayName', 'Curva di Bézier');
    %plot3(x,y,z,'ro', 'DisplayName', 'Punti di interpolazione');
    %legend
    
    figure(2)
    plot3(x_real, y_real, z_real, 'r-', 'DisplayName', 'Funzione interpolata');
    hold on
    plot3(x,y,z,'ro', 'DisplayName', 'Punti di interpolazione');
    plot3(x_interp, y_interp,  z_interp,'b-', 'DisplayName', 'Curva di Bézier');
    plot3(p(1,:),p(2,:),p(3,:),'g+', 'DisplayName', 'Punti di controllo');
    axis([-6 6 -6 6 0 20])
    %legend
end

t_val = linspace(punti(1), punti(end), 3000);
bs_val_err = an_bspl(g, t, t_val);

[x_real_err, y_real_err, z_real_err] = feval(fun, t_val);
x_real_err = x_real_err';
y_real_err = y_real_err';
z_real_err = z_real_err';



tabella = [] %

x_interp_err = bs_val_err*p(1,:)';
y_interp_err = bs_val_err*p(2,:)';
z_interp_err = bs_val_err*p(3,:)';
tabella(k+1, 1) = k;
tabella(k+1, 2) = err;
tabella(k+1, 3) = max_distance([x_interp_err y_interp_err z_interp_err], [x_real_err y_real_err z_real_err]);

while(err>tol && k<100)
    
    k = k+1;
    p = p + w.*d;
    d = [x - (bs*p(1,:)')'; y-(bs*p(2,:)')'; z-(bs*p(3,:)')'];
    
    
    for j=1:length(d)
        norma(j) = norm(d(:,j));
    end
    err = max(norma);
    
    if(k<10 || mod(k, 5)==0)
        fprintf("Errore al passo %1d: %5.5e\n", k, err);
    end
    
    if(grafici > 0)
        [x_real, y_real, z_real] = feval(fun, p_val);
        x_interp = bs_val*p(1,:)';
        y_interp = bs_val*p(2,:)';
        z_interp = bs_val*p(3,:)';
        
        close all
        
        %figure(1)
        %plot3(x_real, y_real, z_real, 'r-', 'DisplayName', 'Funzione interpolata');
        %hold on
        %plot3(x_interp, y_interp,  z_interp,'b-', 'DisplayName', 'Curva di Bézier');
        %plot3(x,y,z,'ro', 'DisplayName', 'Punti di interpolazione');
        %legend
        
        figure(2)
        plot3(x_real, y_real, z_real, 'r-', 'DisplayName', 'Funzione interpolata');
        hold on
        plot3(x,y,z,'ro', 'DisplayName', 'Punti di interpolazione');
        plot3(x_interp, y_interp,  z_interp,'b-', 'DisplayName', 'Curva di Bézier');
        plot3(p(1,:),p(2,:),p(3,:),'g+', 'DisplayName', 'Punti di controllo');
        axis([-6 6 -6 6 0 20])
        %legend
    end
    
    
   
    x_interp_err = bs_val_err*p(1,:)';
    y_interp_err = bs_val_err*p(2,:)';
    z_interp_err = bs_val_err*p(3,:)';
    tabella(k+1, 1) = k;
    tabella(k+1, 2) = err;
    tabella(k+1, 3) = max_distance([x_interp_err y_interp_err z_interp_err], [x_real_err y_real_err z_real_err]);
    
end
legend


%fprintf("Numero iterazioni necessarie: %5.1d\n", k);


if(grafici > 0)
    [x_real, y_real, z_real] = feval(fun, p_val);
    x_interp = bs_val*p(1,:)';
    y_interp = bs_val*p(2,:)';
    z_interp = bs_val*p(3,:)';
    
    close all
    
    %figure(1)
    %plot3(x_real, y_real, z_real, 'r-', 'DisplayName', 'Funzione interpolata');
    %hold on
    %plot3(x_interp, y_interp,  z_interp,'b-', 'DisplayName', 'Curva di Bézier');
    %plot3(x,y,z,'ro', 'DisplayName', 'Punti di interpolazione');
    %legend
    
    figure(2)
    plot3(x_real, y_real, z_real, 'r-', 'DisplayName', 'Funzione interpolata');
    hold on
    plot3(x,y,z,'ro', 'DisplayName', 'Punti di interpolazione');
    plot3(x_interp, y_interp,  z_interp,'b-', 'DisplayName', 'Curva di Bézier');
    plot3(p(1,:),p(2,:),p(3,:),'*', 'DisplayName', 'Punti di controllo');
    %legend
end




