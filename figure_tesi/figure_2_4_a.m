
err = 1;
p = [x; y];
tabella = [];
tic
d = [x - (bs*p(1,:)')'; y-(bs*p(2,:)')'];

%bs_val = an_bspl(g, t, t_real);


for j=1:length(d)
    norma(j) = norm(d(:,j));
end
err = max(norma);
k = 0;

tabella(k+1, 1) = k;
tabella(k+1, 2) = err;
tabella(k+1, 3) = toc;
%fprintf("Errore al passo %1d: %5.5e\n", k, err);

%next_printed = 1;


while(err>tol)
    k = k+1;
    p = p + w.*d;
    d = [x - (bs*p(1,:)')'; y-(bs*p(2,:)')'];
    
    for j=1:length(d)
        norma(j) = norm(d(:,j));
    end
    err = max(norma);
    tabella(k+1, 1) = k;
    tabella(k+1, 2) = err;
    tabella(k+2, 3) = toc;
end
fprintf("Errore al passo %1d: %5.5e\n", k, err);
fprintf("Finito:  %10.10f \n", tabella(k+1, 3));