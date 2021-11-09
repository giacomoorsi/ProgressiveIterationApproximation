% PIA Adattivo
% Lin H. - Adaptive data fitting by the progressive-iterative approximation
%
% parametri da impostare prima della chiamata a pia_bspline_adapt_body:
% x, y                  <- vettori dei punti di interpolazione
% tt                    <- parametri dei punti di interpolazione
% x_real, y_real        <- vettori di punti per disegnare la curva
% t_real                <- parametri dei punti per disegnare la curva
% g                     <- grado della B-Spline
% use_optimal_weigth    <- 1 se WPIA, 0 se PIA
% tol                   <- tolleranza sui punti di interpolazione
% grafici               <- 1 per quello finale, 2 per passo
% adjust_axis           <- 1 per svg
% knot_partition:
%                       1. per prendere i parametri dei punti
%                       2. per punti uniformi in tt(1)...tt(end)
%                       3. per punti  Whitney–Schoenberg




a = min(tt);
b = max(tt);

if (knot_partition == 1)
    t = [repmat(a,1,g), tt, repmat(b,1,g)];
elseif (knot_partition == 2)
    t = [repmat(a,1,g), linspace(a, b, length(tt+2)), repmat(b,1,g)];
elseif (knot_partition == 3)
    t = an_not_a_knot(g,[tt(1) tt tt(end)]);
end
 bs = an_bspl(g, t, tt);

%replico primo ed ultimo punto nel vettore p per avere un numero di punti
%uguale alla dimensione dello spazio spline (n+3)
p = [x(1),x,x(n+1);
    y(1),y,y(n+1)];

% Contiene gli indici dei punti "attivi" cioè sui quali continuo le iterazioni
init_active_points = (2:length(p)-1);
ser_active_points = [];
satisfied_points = [];
unsatisfied_points = [];
%ser_active_points = (2:length(p)-1)
% va da 2 a length(p)-1 perche' il primo e l'ultimo punto sono replicati
% e' necessario traslare di (-1) gli indici se ci si riferisce a posizioni
% nelle matrici d, bs, ecc


if(grafici > 0)
    bs_val = an_bspl(g, t, t_real);
end

if (use_optimal_weight == 1)
    %calcolo autovalori della matrice tridiagonale di iterazione (vedi paper
    %LIN, WANG, DONG, 2003)
    autoval = eig(bs(2:n,3:n+1));
    min_autoval = min(abs(autoval));
    w = 2 / (1 + min_autoval);
else
    w = 1;
end

d = [x - (bs*p(1,:)')'; y-(bs*p(2,:)')'];

for j=1:length(d)
    norma(j) = norm(d(:,j));
end
% disp(norma)
D = max(norma);
k = 0;
fprintf("Errore al passo %2d: %15.6e\n", k, D);


while(length(init_active_points)~=0)
    ser_active_points = init_active_points;
    
    while(length(ser_active_points)~=0)
        
        k = k+1;
        
        % itero solo sui punti ancora attivi
        p(:,ser_active_points) = p(:,ser_active_points) + w.*(d(:,ser_active_points-1));
        
        p(:,1)=p(:,2);
        p(:,n+3)=p(:,n+2);
        
     
        d = [x - (bs*p(1,:)')'; y-(bs*p(2,:)')'];
        
        
        
        % calcolo la norma sugli errori dei punti attivi
        for j=(ser_active_points)
            norma(j-1) = norm(d(:,j-1));
            if(norma(j-1)<tol)
                % rimuovo il punto dai punti attivi
                
                fprintf("Rimosso il punto (%4.5f, %4.5f) per tolleranza sulla norma raggiunta (%4.5e)\n", p(1, j), p(2,j), norma(j-1))
                ser_active_points = ser_active_points(ser_active_points ~= j);
            end
        end
    end
    fprintf("Terminate iterazioni. Procedo alla classificazione. \n");
    
    % compute the difference vector
    
    indices = union(init_active_points, satisfied_points);
    indices = union(indices, unsatisfied_points);
    if (size(indices, 1) > size(indices, 2))
        indices = indices';
    end
    
    norma2 = [];
    init_active_points = [];
    eps_0 = tol;
    d_g = [];
    d_l = [];
    d_l_ser = [];
    i=1;
    for i=indices
        norma2(i) = norm(d(:,i-1));
        if norma2(i) > eps_0
            d_g(:, end+1) = d(:, i-1);
            init_active_points(end+1) = i;
        else
            d_l(:,end+1) = d(:, i-1);
            d_l_ser(end+1) = i;
        end
    end
    
    D = max(norma2);
    fprintf("Errore di interpolazione al passo %2d: %15.6e\n", k, D);

    % categorize d_l into subset S_p where each difference vector satisfies
    % Eq. 5 and subset U_p where each vector does not satifly Eq. 5 
    
    fprintf("Punti ancora attivi: ");
    disp(init_active_points);
    
    if (length(init_active_points)~=0)
        b_k = bs(init_active_points-1, init_active_points-1);
    
        precalc = 1/(1-norm(eye(length(b_k)) - b_k,Inf));
        
        mnorm = norm(norma2(init_active_points),Inf);
        
        s_p = [];
        u_p = [];
        
        [n1,n2]=size(d_l);        
        for (i=1:n2)
            L = bs(init_active_points-1, d_l_ser(i))';
            if (d_l(i) <= eps_0 - precalc * norm(L,Inf) * mnorm)
                s_p(end+1) = d_l_ser(i);
            else
                u_p(end+1) = d_l_ser(i);
            end
        end
        
        if length(u_p)~=0
            satisfied_points = s_p;
            unsatisfied_points = u_p;
        else
            unsatisfied_points = [];
            satisfied_points = [];
        end
       
        disp('unsatisfied points')
        disp(unsatisfied_points)
        
    else
        fprintf("init_active_points vuoto! \nAlgoritmo terminato! \n");
    end
    
end

fprintf("Errore di interpolazione al passo %2d: %15.6e\n", k, D);

if(grafici > 0)
    x_interp = bs_val*p(1,:)';
    y_interp = bs_val*p(2,:)';
    
    figure(2)
    plot(x_real, y_real, 'r-', 'DisplayName', 'Curva interpolata');
    hold on
    
    if (adjust_axis)
        axis ij
        axis equal    
    end
    
    plot(x_interp, y_interp, 'b-', 'DisplayName', 'B-Spline');
    plot(x,y,'ro', 'DisplayName', 'Punti di interpolazione');
    legend;
    
    figure(3)
    plot(x_real, y_real, 'r-', 'DisplayName', 'Curva interpolata');
    hold on
    
    if (adjust_axis)
        axis ij
        axis equal    
    end
    plot(x_interp, y_interp, 'b-', 'DisplayName', 'B-Spline');
    plot(p(1,:),p(2,:),'*', 'DisplayName', 'Punti di controllo');
    legend;
end