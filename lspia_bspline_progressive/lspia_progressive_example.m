clear all
g = 3;
grafici = 2;

n_q = 100;
n = 10;
knot_partition = 3;

R = loadsvg('clef_lspia.svg', 0.01,0 );
C = R{1,1}';
C = rescale(C);
parametri = curv2_param(1, C(1,:), C(2,:));

equidistanti = 1;

% % campiono n_q punti dall'immagine SVG
if(equidistanti == 0)
    punti = []
    i=0;
    for k = floor(linspace(1, length(C), n_q))
        punti(i+1) = k;
        q_x(i+1) =C(1,k);
        q_y(i+1) =C(2,k);
        i = i+1;
    end
    q_tt = curv2_param(1, q_x, q_y);
    
    a = min(punti);
    b = max(punti);
    
else
    [punti, indici] = curvspace(C', n_q);
    q_x = punti(:,1)';
    q_y = punti(:,2)';
    q_tt = curv2_param(1, q_x, q_y);
    a = min(indici);
    b = max(indici);
   % punti = floor(linspace(a, b, n_q));
   %punti = indici;
  
  
end



% prendo n+1 punti tra quelli campionati
indici = floor(linspace(1, length(punti), n+1));


x = q_x(indici);
y = q_y(indici);
tt = q_tt(indici);

tol1 = 1e-5; % tolleranza nei punti
tol2 = 1e-9; % avanzamento minimo

adjust_axis = 1;


iterations_before_knot_insertion = 30; % numero di iterazioni da eseguire prima di inserire un nuovo nodo
max_knots_inserted = 1000; % numero massimo di nodi da inserire



max_iterations = 10000; % numero massimo di iterazioni



lspia_progressive_bspline_body;

