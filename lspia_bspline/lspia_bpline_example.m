clear all
clc

R = loadsvg('clef_lspia.svg', 0.01,0 );
C = R{1,1}';
parametri = curv2_param(1, C(1,:), C(2,:)); % calcolo parametri dei punti campionati

g = 3; % grado della spline
grafici = 3; % 3 per mostrare tutti i grafici
n_q = 70; % numero punti di approssimazione
n = 15; % numero punti di controllo
equidistanti = 1; % 1 per estrarre punti equidistanti dai punti campionati

knot_partition = 3; % 1 per equistanziati, 2 per parametri, 3 per Whitney–Schoenberg



if(equidistanti == 0)
    punti = []
    i=0;
    for k = floor(linspace(1, length(C), n_q))
        punti(i+1) = k;
        q_x(i+1) =C(1,k);
        q_y(i+1) =C(2,k);
        i = i+1;
    end
    
    a = min(punti);
    b = max(punti);
    
else
    [punti, indici] = curvspace(C', n_q);
    q_x = punti(:,1)';
    q_y = punti(:,2)';
    q_tt = parametri(indici); % tengo solo i parametri dei punti utilizzati
    
end

% prendo n+1 punti tra quelli campionati (saranno i punti iniziali LSPIA)
indici = floor(linspace(1, n_q, n+1));

x = q_x(indici);
y = q_y(indici);
tt = q_tt(indici);

% tt : parametri dei punti x;y
% punti: parametri dei punti 


tol1 = 1e-5; % tolleranza nei punti
tol2 = 1e-11; % avanzamento minimo

adjust_axis = 1; % 1 per gli svg

lspia_bspline_body;

