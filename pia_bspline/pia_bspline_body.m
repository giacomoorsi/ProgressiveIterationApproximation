% PIA / WPIA su B-Spline
%
% parametri da impostare prima della chiamata a pia_bspline_body:
% x, y                  <- vettori dei punti di interpolazione
% tt                    <- parametri dei punti di interpolazione
% x_real, y_real        <- vettori di punti per disegnare la curva
% t_real                <- parametri dei punti per disegnare la curva
% g                     <- grado della B-Spline
% use_optimal_weigth    <- 1 se WPIA, 0 se PIA
% tol                   <- tolleranza sui punti di interpolazione
% grafici               <- 1 per mostrarli
% adjust_axis           <- 1 per svg
% knot_partition:
%                       1. per prendere i parametri dei punti
%                       2. per punti uniformi in tt(1)...tt(end)
%                       3. per punti  Whitney–Schoenberg



p = [x;y];
a = min(tt);
b = max(tt);


k = length(p)-g-1;
if (knot_partition == 1)
    t = [repmat(a,1,g+1), tt(3:end-2), repmat(b,1,g+1)];
elseif (knot_partition == 2)
    t = [repmat(a,1,g), linspace(a, b, k+2), repmat(b,1,g)];
elseif (knot_partition == 3)
    t = an_not_a_knot(g,tt);
end
bs = an_bspl(g, t, tt);


if (use_optimal_weight == 1)
    autoval = eig(bs);
    min_autoval = min(abs(autoval));
    w = 2 / (1 + min_autoval);
else
    w = 1;
end

d = [x - (bs*p(1,:)')'; y-(bs*p(2,:)')'];

bs_val = an_bspl(g, t, t_real);


for j=1:length(d)
    norma(j) = norm(d(:,j));
end
err = max(norma);
k = 0;
fprintf("Errore al passo %1d: %5.5e\n", k, err);

next_printed = 1;


while(err>tol)
    k = k+1;
    p = p + w.*d;
    d = [x - (bs*p(1,:)')'; y-(bs*p(2,:)')'];
    
    for j=1:length(d)
        norma(j) = norm(d(:,j));
    end
    err = max(norma);
    
    if(k==next_printed)
        next_printed = next_printed * 2;
        fprintf("Errore al passo %1d: %5.5e\n", k, err);
    end
end
fprintf("Errore al passo %1d: %5.5e\n", k, err);

if(grafici > 0)
    x_interp = bs_val*p(1,:)';
    y_interp = bs_val*p(2,:)';
    
    figure(1)
    plot(x_real, y_real, 'r-', 'DisplayName', 'Curva interpolata');
    hold on
    plot(x_interp, y_interp, 'b-', 'DisplayName', 'B-Spline');
    plot(x,y,'ro', 'DisplayName', 'Punti di interpolazione');
    legend
    if (adjust_axis)
       axis ij
       axis equal
    end
        
    
    figure(2)
    plot(x_real, y_real, 'r-', 'DisplayName', 'Curva interpolata');
    hold on
    plot(x_interp, y_interp, 'b-', 'DisplayName', 'B-Spline');
    plot(p(1,:),p(2,:),'*', 'DisplayName', 'Punti di controllo');
    legend
    if (adjust_axis)
       axis ij
       axis equal
    end
end