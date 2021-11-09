% PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = nrbline ([0 0], [1 0]);

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2];

% Physical parameters
problem_data.c_diff  = @(x) ones(size(x));
problem_data.grad_c_diff = @(x) reshape (zeros(size(x)), [1, size(x)]);

% Source and boundary terms
problem_data.f = @(x) (2*pi)^2*sin(2*pi*x);
problem_data.h = @(x, ind) zeros (size (x));

% Exact solution (optional)
problem_data.uex     = @(x) sin(2*pi*x);
problem_data.graduex = @(x) 2*pi*cos(2*pi*x);

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = 2;            % Degree of the splines
method_data.regularity  = 1;            % Regularity of the splines
method_data.nsub_coarse = 2;            % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = 2;            % Number of subdivisions for each refinement
method_data.nquad       = 3;            % Points for the Gaussian quadrature rule
method_data.space_type  = 'simplified'; % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 0;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
% adaptivity_data.flag = 'elements';
adaptivity_data.flag = 'functions';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = .5;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 12;
adaptivity_data.max_ndof = 5000;
adaptivity_data.num_max_iter = 15;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-9;

% GRAPHICS
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

[geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);

% Plot
[eu, F] = sp_eval (u, hspace, geometry, 101);
figure; plot (F, eu)


%!test
%! problem_data.geo_name = nrbline ([0 0], [1 0]);
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2];
%! problem_data.c_diff  = @(x) ones(size(x));
%! problem_data.grad_c_diff = @(x) reshape (zeros(size(x)), [1, size(x)]);
%! problem_data.f = @(x) (2*pi)^2*sin(2*pi*x);
%! problem_data.h = @(x, ind) zeros (size (x));
%! problem_data.uex     = @(x) sin(2*pi*x);
%! problem_data.graduex = @(x) 2*pi*cos(2*pi*x);
%! method_data.degree      = 2;            % Degree of the splines
%! method_data.regularity  = 1;            % Regularity of the splines
%! method_data.nsub_coarse = 2;            % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
%! method_data.nsub_refine = 2;            % Number of subdivisions for each refinement
%! method_data.nquad       = 3;            % Points for the Gaussian quadrature rule
%! method_data.space_type  = 'simplified'; % 'simplified' (only children functions) or 'standard' (full basis)
%! method_data.truncated   = 0;            % 0: False, 1: True
%! adaptivity_data.flag = 'functions';
%! adaptivity_data.C0_est = 1.0;
%! adaptivity_data.mark_param = .5;
%! adaptivity_data.mark_strategy = 'MS';
%! adaptivity_data.max_level = 12;
%! adaptivity_data.max_ndof = 5000;
%! adaptivity_data.num_max_iter = 11;
%! adaptivity_data.max_nel = 5000;
%! adaptivity_data.tol = 1e-9;
%! plot_data.print_info = false;
%! plot_data.plot_hmesh = false;
%! plot_data.plot_discrete_sol = false;
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.iter, 11)
%! assert (solution_data.ndof, [4 6 10 18 34 58 98 122 218 246 430]);
%! assert (solution_data.nel, [2 4 8 16 32 56 96 120 216 244 428]);
%! assert (solution_data.err_h1s, [6.46392167473527e-02, 5.35889596011715e-01, 1.09845482621112e-01, 2.60112950295797e-02,...
%!          6.41325455175922e-03, 1.97502613294894e-03, 7.98167327328958e-04, 4.12338369102373e-04, 1.32998923347470e-04, 1.01034043201185e-04, 3.31687093456353e-05], 1e-15)
%! assert (solution_data.flag, 2)
%!
%! adaptivity_data.flag = 'elements';
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.iter, 11)
%! assert (solution_data.ndof, [4 6 10 12 22 30 54 62 110 122 210])
%! assert (solution_data.nel, [2 4 8 12 20 28 52 60 108 120 208]);
%! assert (solution_data.err_h1s, [6.46392167473527e-02   5.35889596011715e-01   1.09845482621112e-01   8.94046382886190e-02   2.09289645926629e-02   8.31121898502044e-03   2.55125993516858e-03   1.66523666575394e-03 5.41035873522990e-04   4.12338369102373e-04   1.47630836012780e-04], 1e-15)
%! assert (solution_data.flag, 2)
%!
%! method_data.space_type  = 'standard';
%! adaptivity_data.flag = 'functions';
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.iter, 11)
%! assert (solution_data.ndof, [4 6 10 18 34 58 98 122 218 246 430]);
%! assert (solution_data.nel, [2 4 8 16 32 56 96 120 216 244 428]);
%! assert (solution_data.err_h1s, [6.46392167473527e-02, 5.35889596011715e-01, 1.09845482621112e-01, 2.60112950295797e-02,...
%!          6.41325455175922e-03, 1.97502613294894e-03, 7.98167327328958e-04, 4.12338369102373e-04, 1.32998923347470e-04, 1.01034043201185e-04, 3.31687093456353e-05], 1e-15)
%! assert (solution_data.flag, 2)
%!
%! method_data.space_type  = 'standard';
%! adaptivity_data.flag = 'elements';
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.iter, 11)
%! assert (solution_data.ndof, [4 6 10 14 18 30 54 62 110 122 210])
%! assert (solution_data.nel, [2 4 8 12 16 28 52 60 108 120 208]);
%! assert (solution_data.err_h1s, [6.46392167473527e-02   5.35889596011715e-01   1.09845482621112e-01   6.42010867661501e-02   2.60112950295797e-02   8.31121898502044e-03   2.55125993516858e-03   1.66523666575394e-03 5.41035873522990e-04   4.12338369102373e-04   1.47630836012780e-04], 1e-15)
%! assert (solution_data.flag, 2)
%!
%! adaptivity_data.max_level = 4;
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.flag, 3)
%! assert (solution_data.iter, 4)
%!
%! adaptivity_data.max_level = 12;
%! adaptivity_data.max_ndof = 50;
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.flag, 4)
%! assert (solution_data.iter, 7)
%!
%! adaptivity_data.max_ndof = 5000;
%! adaptivity_data.max_nel = 50;
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.flag, 5)
%! assert (solution_data.iter, 7)
%!
%! adaptivity_data.flag = 'functions';
%! adaptivity_data.max_nel = 5000;
%! adaptivity_data.tol = 1e-2;
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.flag, 1)
%! assert (solution_data.iter, 7)
