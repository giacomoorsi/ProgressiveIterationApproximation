close all
clear all


R = loadsvg('clef_lspia.svg', 0.01,0 );
C = R{1,1}';

% calcolo la parametrizzazione dei punti
parametri = curv2_param(1, C(1,:), C(2,:));

g = 3; % grado B-Spline
grafici = 1; % grafici
n = 8;
equidistanti = 0; % 1 per estrarre punti equidistanti dai punti campionati
use_optimal_weight = 0;
tol = 1e-10; 
knot_partition = 3; % 1 per equistanziati, 2 per parametri, 3 per Whitney–Schoenberg
adjust_axis = 1;

p_real = floor(linspace(1, length(C), 200));
t_real = parametri(p_real);
x_real = C(1,p_real);
y_real = C(2,p_real);

if(equidistanti == 0)
    punti = []
    i=0;
    for k = floor(linspace(1, length(C), n))
        punti(i+1) = k;
        x(i+1) =C(1,k);
        y(i+1) =C(2,k);
        i = i+1;
    end
    
    tt = parametri(punti);
    
else
    [punti, indici] = curvspace(C', n);
    x = punti(:,1)';
    y = punti(:,2)';
    tt = parametri(indici); % tengo solo i parametri dei punti utilizzati
end

adjust_axis = 1; % 1 per gli svg

pia_bspline_body;

