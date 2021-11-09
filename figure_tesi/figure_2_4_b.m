tic
p_x = bs\x';
p_y = bs\y';
fine = toc;
fprintf("finito %10.10f", fine)


%t_val = linspace(a, b, 500);
%bs_val = an_bspl(g, t, t_val);

hold off
plot(tabella(:, 3), tabella(:, 2), 'k-');
hold on
plot([0; fine], [1e-1; 1e-1], 'g-')
plot([fine; fine], [1e-1; endd], 'g-')
plot([fine; end_time], [endd; endd], 'g-')
set(gca, 'YScale', 'log')