%Geometry domain
L=20; H=1; B=1; A=H*B; I = B*H^3/12;
interval = nrbline ([0 0],[L 0]);
problem_data.geo_name = interval;

%Boundary conditions
%Left end: Moment + concentrated force
%Right end: clamping
problem_data.drchlt1_ends = [false true]; %deflections, homogeneous
problem_data.drchlt2_ends = [false true]; %rotation angles, homogeneous
problem_data.nmnn1_ends   = [true false]; %forces 
problem_data.nmnn2_ends   = [true false]; %moments

M(1) = -20; M(2) = 0;
P = 2;
problem_data.M = M; %moment values
problem_data.P = P; %force value 

f_c = 0.5;
problem_data.f = @(x) f_c * sin(pi*x/L); %distributed loading

%Physical parameters
E   = 210000;
problem_data.EI = @(x) E * I* ones (size (x));

% Discretization parameters
method_data.degree     = 3;     % Degree of the basis functions
method_data.regularity = 2;     % Regularity of the basis functions
method_data.nsub       = 20;     % Number of subdivisions
method_data.nquad      = 5;     % Points for the Gaussian quadrature rule

% Solving
[geometry, msh, space, w] = solve_BE_Beam (problem_data, method_data);

%Postprocessing
%Exact solution
c1 = (L*f_c+P*pi)/(pi*E*I);
c2 = M(1)/(E*I);
c3 = -L*(L^2*pi^2*f_c+L*P*pi^3+2*M(1)*pi^3-2*L^2*f_c)/(pi^3*E*I*2);
c4 = L^2*(2*L^2*pi^2*f_c+2*L*P*pi^3+3*M(1)*pi^3-6*L^2*f_c)/(pi^3*E*I*6);

w_ex = @(x) f_c*sin(pi*x/L)/E/I*(L/pi)^4 + c1*x.^3/6 + c2*x.^2/2 + c3*x + c4;
%Gradient of exact solution
gradwex = @(x) f_c*cos(pi*x/L)/E/I*(L/pi)^3 + c1*x.^2/2 + c2*x + c3;
%Hessian of exact solution
hesswex = @(x) -f_c*sin(pi*x/L)/E/I*(L/pi)^2 + c1*x + c2;

%L2, H1 and H2 norms of error                 

[errh2, errh1, errl2] = sp_h2_error (space, msh, w, w_ex, gradwex, hesswex);
errors = [errh2, errh1, errl2]
[normh2, normh1, norml2] = sp_h2_error (space, msh, zeros (size (w)), w_ex, gradwex, hesswex);
relative_errors = [errh2/normh2, errh1/normh1, errl2/norml2]

[eu, F] = sp_eval (w, space, geometry, 100);
han = figure;
plot(F,eu,F,w_ex(F),'--')
legend('FEM-solution', 'Exact solution')
xlabel('x')
ylabel('w')

%!test
%! ex_BE_beam_static_2
%! close (han)
%! assert (errors, [3.42006802485252e-06 5.44598574884688e-07 1.51095034097015e-07], 1e-11);
%! assert (relative_errors, [3.46016655651430e-06 5.51007520112053e-07 1.53624073269057e-07], 1e-11);
