close all
clear all
clc

knot_partition = 3;
n = 20;
g = 3;
fun = 'gerono';
use_optimal_weight = 0;
tol = 1e-15;
grafici = 1;
adjust_axis = 0;

a = 0;
b = 2*pi;
tt = linspace(0, 2*pi, n+1);
[x,y] = feval(fun, tt);

t_real = linspace(tt(1), tt(end), 100);
[x_real, y_real] = feval(fun, t_real);

pia_bspline_adapt_body;