
clear all
g = 3;
grafici = 2;

n_q = 300;
n = 29;
knot_partition = 3;

R = loadsvg('rapace.svg', 0.2,0 );
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

q_plot = linspace(q_tt(1), q_tt(end), 1000);
    


% prendo n+1 punti tra quelli campionati
indici = floor(linspace(1, length(punti), n+1));


x = q_x(indici);
y = q_y(indici);
tt = q_tt(indici);

tol1 = 5e-8; % tolleranza nei punti
tol2 = 1e-5; % avanzamento minimo

adjust_axis = 0;


iterations_before_knot_insertion = 30; % numero di iterazioni da eseguire prima di inserire un nuovo nodo
max_knots_inserted = 700; % numero massimo di nodi da inserire



max_iterations = 10000; % numero massimo di iterazioni

% LSPIA progressive su B-Spline
%
%
% parametri da impostare prima della chiamata a lspia_progressive_bspline_body:
% q_x, q_y              <- vettori punti di approssimazione
% q_tt                  <- parametri dei punti di approssimazione
% x, y                  <- vettori punti di approssimazione iniziali
% tt                    <- parametri dei punti x, y
% tol1                  <- tolleranza sull'errore di approssimazione
% tol2                  <- tolleranza sull'avanzamento (diff. errori al
%                       passo k e al passo k-1)
% grafici               <- 1 per grafico conclusivo, 2 per passo passo
% knot_partition:
%                       1. per prendere i parametri dei punti
%                       2. per punti uniformi in tt(1)...tt(end)
%                       3. per punti  Whitney???Schoenberg
% g                     <- grado della B-Spline
% adjust_axis           <- 1 per svg

% Variabili facoltative:
% iterations_before_knot_insertion      <- numero di iterazioni da eseguire
%                                       prima di inserire un nuovo nodo
%                                       (default 50)
% max_knots_inserted                    <- numero massimo di nodi da inserire
%                                       (default 30)
% max_iterations                        <- numero massimo di iterazioni
%                                       (default 500)
% fun                                   <- contiene il nome della funzione
%                                       che ?? possibile valutare per mostrare la
%                                       funzione valutata nei nuovi nodi



% L'algoritmo itera "iterations_before_knot_insertion" volte prima di
% aggiungere un nuovo nodo e si interrompe quando si verifica una delle
% seguenti condizioni:
% a) ?? stato superato il numero massimo di iterazioni (max_iterations)
% b) l'avanzamento (differenza tra l'errore di approssimazione del passo
% attuale con il precendente) ?? minore di tol2
% c) l'errore di approssimazione massimo ?? minore di tol1
% d) sono stati inseriti max_knots_inserted nodi e completate le
% iterations_before_knot_insertion iterazioni successive


% numero di iterazioni da eseguire prima di inserire un nuovo nodo
if(~exist('iterations_before_knot_insertion'))
    iterations_before_knot_insertion = 50;
end

% numero massimo di nodi da inserire
if(~exist('max_knots_inserted'))
    max_knots_inserted = 30;
end

% numero massimo di iterazioni
if(~exist('max_iterations'))
    max_iterations = 500;
end


p = [x;y];


k = length(p)-g-1;
if (knot_partition == 1)
    t = [repmat(a,1,g+1), tt(3:end-2), repmat(b,1,g+1)];
elseif (knot_partition == 2)
    t = [repmat(a,1,g), linspace(a, b, k+2), repmat(b,1,g)];
elseif (knot_partition == 3)
    t = an_not_a_knot(g,tt);
end
bs = an_bspl(g, t, q_tt);
bs_plot = an_bspl(g, t, q_plot);

btb = bs'* bs;

% Peso esatto
% autoval = eig(btb);
% mu = 2/(min(autoval)+max(autoval));

% Peso ottimizzato
C = max(sum(btb, 2)); % prendo il massimo della somma di ogni riga
mu = 2/C; % definisco mu come da approssimazione


err1 = tol1 + 1;
avanzamento = tol2 + 1;
h = 1; % contatore iterazioni dall'ultimo knot insertion
knots_inserted = 0; % numero nodi gi?? inseriti
k = 1; % contatore iterazioni totali
err_interp_max(k) = tol1+1;
err_approx_max = tol2+1;

indice_p_nodo = zeros(1, length(t)); % nella posizione i-esima contiene l'indice del primo punto appartenente all'intervallo nodale [t[i], t[i+1])
j=1;
i=1;
while(i<=length(q_tt) && j<=length(t))
    while(q_tt(i)>=t(j+1))
        j = j+1;
        if(j==length(t))
            break
        end
    end
    indice_p_nodo(j) = i;
    i = i+1;
    while(i<=length(q_tt) && q_tt(i)<t(j+1))
        i=i+1;
    end
end
indice_p_nodo(length(t)-g) = indice_p_nodo(end);
indice_p_nodo(end) = 0;

tabella = [];

h = 0;


% && avanzamento>=tol2
while(k<max_iterations  && err_approx_max(k)>=tol1 && ~(h==iterations_before_knot_insertion-1 && knots_inserted==max_knots_inserted))
    h = h+1;
    k = k+1;
    
    x_approx = (bs*p(1,:)')';
    y_approx = (bs*p(2,:)')';
    delta = [q_x - x_approx; q_y - y_approx];
    for(i=1:(n+1+knots_inserted))
        adj(1, i) = mu * sum(bs(:, i)'.*delta(1, :));
        adj(2, i) = mu * sum(bs(:, i)'.*delta(2, :));
    end
    p = p + adj;
    
    
    for j=1:length(delta)
        norma(j) = norm(delta(:,j));
    end
    err_approx_max(k) = max(norma);
    tabella(k,1) = k-2;
    tabella(k,2) = length(p);
    tabella(k,3) = err_approx_max(k);
    
    avanzamento = abs(err1 - err_approx_max(k));
    err1 = err_approx_max(k);
    
%     if((k<10 || mod(k, 5)==0) && grafici > 0)
%         fprintf("Errore al passo %1d: %5.5e\n", k, err1);
%         fprintf("Avanzamento al passo %1d: %5.5f\n", k, avanzamento);
%     end
    
    
    
    if((k<10 || mod(h, 5)==0)  && grafici == 2)
        a = 1;
            figure(1)
            hold off
            plot(q_x, q_y, 'b.', 'DisplayName', 'Punti di approssimazione');
            hold on
           % axis ij
            axis equal
            x_plot = (bs_plot*p(1,:)')';
            y_plot = (bs_plot*p(2,:)')';

            plot(x_plot, y_plot, 'r-', 'DisplayName', strcat('B-Spline k:', num2str(k)));
            plot(p(1,:), p(2,:),'g+:', 'DisplayName', 'Punti di controllo');
             %legend;
            %axis([-0.4 0.85 -0.25 1.25])
            axis([-0.1 0.9 -0.1 1.1])
            %title(strcat('k: ', num2str(k), ", h:", num2str(h)));
            fprintf("k: %1.0f, numero nodi: %1.0f, h: %1.0f \n", k, length(p), h)
        end

    
    % inserisco un nuovo nodo dopo aver fatto
    % iterations_before_knot_insertion senza inserire nodi
    if(h == iterations_before_knot_insertion)
        h = 0;
        
        
        
        
        
        % mostro segmento
%         if (grafici == 2)
%             hold on
%             for i=1:length(q_x)
%                 plot([x_approx(i), q_x(i)],[y_approx(i), q_y(i)])
%             end
%         end
        
        % calcolo per ogni partizione nodale la quantit?? di errori in essa
        %contenuti
        err_approx_interval = zeros(1, length(t));
        j=1;
        i=1;
        while(j<=length(t) && i<=length(q_tt))
            while(j<length(t) && q_tt(i)>=t(j+1))
                j = j+1;
            end
            err_approx_interval(j) = err_approx_interval(j) + norma(i);
            i = i+1;
        end
        
        % calcolo l'intervallo nel quale si verifica l'errore massimo
        [err_approx_interval_max, new_knot_index] = max(err_approx_interval);
        
        
        knots_inserted = knots_inserted + 1;
        
        % per calcolare il newknot uso la formula di  LSPIA progressive
        % (cap. 5.2 LSPIA)
        err_parziale = 0;
        
        % cerco il parametro del primo punto che appartiene all'intervallo
        % BACKUP: vecchia versione
%                 i=1;
%                 while(q_tt(i)<t(new_knot_index))
%                    i=i+1;
%                 end
%                 fprintf("Trovato %5.4f", i);
        
        i=indice_p_nodo(new_knot_index);
        
        % scorro i punti dell'intervallo finch?? l'errore ?? minore di
        % err_approx_interval_max/2. Quando mi fermo ho trovato il punto
        % esatto nel quale conviene aggiungere un nodo
        while(err_parziale<err_approx_interval_max/2 && i<length(q_tt)-1)
            err_parziale = err_parziale + norma(i);
            i = i+1;
        end
        
        if (i == length(q_tt))
            newknot = (q_tt(end-1) + q_tt(end))/2;
        else
            newknot = (q_tt(i) + q_tt(i+1))/2;
        end
        
        
        % knot insertion
        [p, t] =  bspkntins(g,p,t,newknot);
        bs = an_bspl(g,t,q_tt); % ottimizzabile
        if(grafici == 2)
            
           bs_plot = an_bspl(g,t,q_plot);
        end
        fprintf("Aggiunto nuovo nodo in posizione %1.0d con valore %4.5f\n", new_knot_index, newknot)
        fprintf("Passo %2.0f, errore: %5.5e\n", k, err_approx_max(k))
        % aggiorno indice_p_nodo
        % fino alla posizione new_knot_index deve rimanere tutto inviato.
        % Nella posizione new_knot_index+1 inserisco i+1 ovvero il primo
        % nodo del nuovo intervallo. Rimane da aggiornare soltanto la
        % posizione new_knot_index+2, per la quale bisogna cercare quale
        % sia il primo punto dell'intervallo
        indice_p_nodo = [indice_p_nodo(1:new_knot_index) (i+1) indice_p_nodo(new_knot_index+1:end)];
        if(new_knot_index+2 < length(indice_p_nodo)-g)
            while(q_tt(i) < t(new_knot_index+2))
                i = i+1;
            end
            indice_p_nodo(new_knot_index+2) = i;
        end
        
        
        % Codice di backup che ricalcola tutto indice_p_nodo scorrendo tutti i punti. 
        % E' possibile vedere che i risultati che danno sono identici
%         indice_p_nodo_back = zeros(1, length(t)); % nella posizione i-esima contiene l'indice del primo punto appartenente all'intervallo nodale [t[i], t[i+1])
%         j=1;
%         i=1;
%         while(i<=length(q_tt) && j<=length(t))
%             while(q_tt(i)>=t(j+1))
%                 j = j+1;
%                 if(j==length(t))
%                     break
%                 end
%             end
%             indice_p_nodo_back(j) = i;
%             i = i+1;
%             while(i<=length(q_tt) && q_tt(i)<t(j+1))
%                 i=i+1;
%             end
%         end
%         
%         indice_p_nodo_back(length(t)-g) = indice_p_nodo_back(end);
%         indice_p_nodo_back(end) = 0;
        

        btb = bs'* bs;

        % ricalcolo peso
        C = max(sum(btb, 2));
        mu = 2/C;
        
        
        
        if (exist('fun') && grafici == 2)
            % valuto la funzione nel nuovo punto per mostrarlo sul grafico
            
            [newp_val_x, newp_val_y] = feval(fun, newknot);
            figure(1)
            plot(newp_val_x, newp_val_y, '*', 'DisplayName', 'Funzione valutata nel nuovo nodo')
            pause
            fprintf("Aggiunto nuovo nodo (%4.5f) in posizione (%4.0f). Funzione valutata nel punto: (%4.5f, %4.5f)  \n", newknot, new_knot_index,newp_val_x, newp_val_y)
            
        end
        
       
    end
end

if (k>=max_iterations)
   fprintf("Interrotto per massime iterazioni raggiunte: %1.0f \n", k); 
end

if ( avanzamento<tol2)
   fprintf("Interrotto per avanzamento: %5.5e \n",avanzamento); 
end

if ( err_approx_max(k)<tol1)
   fprintf("Interrotto per tolleranza errore approssimazione: %5.5e \n",err_approx_max(k)); 
end

if ((h==iterations_before_knot_insertion-1 && knots_inserted==max_knots_inserted))
    fprintf("Interrotto numero max nodi inseriti : %5.0f \n",knots_inserted); 
end

if (grafici > 0)
    figure(1)
            hold off
            plot(q_x, q_y, 'b.', 'DisplayName', 'Punti di approssimazione');
            hold on
           % axis ij
            axis equal
             bs_plot = an_bspl(g,t,q_plot);
            x_plot = (bs_plot*p(1,:)')';
            y_plot = (bs_plot*p(2,:)')';

            plot(x_plot, y_plot, 'r-', 'DisplayName', strcat('B-Spline k:', num2str(k)));
            plot(p(1,:), p(2,:),'g+:', 'DisplayName', 'Punti di controllo');
             %legend;
            %axis([-0.4 0.85 -0.25 1.25])

            %title(strcat('k: ', num2str(k), ", h:", num2str(h)));
            fprintf("k: %1.0f, numero nodi: %1.0f, h: %1.0f \n", k, length(p), h)
end