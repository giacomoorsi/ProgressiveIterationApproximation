
clear all
clc

knot_partition = 3;
n = 2000;
g = 3;
fun = 'epitrochoid';
use_optimal_weight = 0;
tol = 1e-15;
grafici = 0;
adjust_axis = 0;


a = 0;
b = 6*pi;
tt = linspace(a, b, n);
[x, y] = feval(fun, tt);
fprintf("punti calcolati \n");


%t_real = linspace(tt(1), tt(end), 500);
%[x_real, y_real] = feval(fun, t_real);

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
fprintf("parto \n");



function [x,y] = epitrochoid(t)
% definita tra 0 e 6*pi
a=2;
b=3;
c=1;
x = (a-b)*cos(t)-c*cos((a/b+1)*t);
y = (a-b)*sin(t)-c*sin((a/b+1)*t);
end