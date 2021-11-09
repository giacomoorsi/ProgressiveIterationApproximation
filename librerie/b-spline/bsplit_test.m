
% definisco 

% 3 curve di Bézier cubiche
k = [1,2]; % nodi

n = 2;

a = 0;
b = 1;

% il numero di punti di controllo è n+k+1 = 6


% punti di controllo
% P = [0, 0;
%     0.5, 7; 
%     1, 3;
%     1.5, 10;
%     2, 0; 
%     2.5, 5  ;
%     ];
P = [0.25, 0.1;
    0.1, 0.4;
    0.3, 0.9;
    0.6, 0.7;
    0.5, 0.4;
    0.7, 0.1;
    0.9, 0.3;]


% partizione nodale
k = [0.20, 0.40, 0.60, 0.8];


%t = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3];
t = [0, 0, 0, 0.20, 0.40, 0.6, 0.8, 1, 1, 1];

length(t)
% length(x) = 2n+k+2 = 6+2+2 = 10

%dimensione spazio spline
mpk=length(t)-n-1

%for x=0:0.01:3
%    y = an_bspl(n,t,x)*P;
    
%end

%x = linspace(0, 3, 100);

x = linspace(a, b, 100);
Y = an_bspl(n,t,x)*P;
figure(1)
plot(Y(:, 1), Y(:, 2));
hold on
plot(P(:, 1), P(:, 2), '*-');









