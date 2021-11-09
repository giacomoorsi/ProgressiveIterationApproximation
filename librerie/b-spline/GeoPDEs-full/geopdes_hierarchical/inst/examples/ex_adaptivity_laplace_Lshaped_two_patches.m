% PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_Lshaped_two_patches.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));
problem_data.grad_c_diff = @(x, y) cat (1, ...
            reshape (zeros(size(x)), [1, size(x)]), ...
            reshape (zeros(size(x)), [1, size(x)]));

% Singular function
k = 1; % Constant that characterizes the singularity
problem_data.f = @(x, y) zeros (size (x));
problem_data.h = @(x, y, ind) singular_function_laplace (x, y, k);

% Exact solution (optional)
problem_data.uex     = @(x, y) singular_function_laplace (x, y, k);
problem_data.graduex = @(x, y) singular_function_maxwell (x, y, k);

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [3 3];        % Degree of the splines
method_data.regularity  = [2 2];        % Regularity of the splines
method_data.nsub_coarse = [2 2];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.nquad       = [4 4];        % Points for the Gaussian quadrature rule
method_data.space_type  = 'simplified'; % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 0;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
% adaptivity_data.flag = 'elements';
adaptivity_data.flag = 'functions';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = .5;
adaptivity_data.mark_strategy = 'MS';
adaptivity_data.max_level = 10;
adaptivity_data.max_ndof = 15000;
adaptivity_data.num_max_iter = 5;
adaptivity_data.max_nel = 15000;
adaptivity_data.tol = 1e-10;

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

[geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);

% EXPORT VTK FILE
npts = [51 51];
output_file = 'laplace_adaptivity_Lshaped.vts';
sp_to_vtk (u, hspace, geometry, npts, output_file, {'solution', 'gradient', 'laplacian'}, {'value', 'gradient', 'laplacian'})

plot_numerical_and_exact_solution(u,hspace,geometry,npts,problem_data.uex)

%!test
%! problem_data.geo_name = 'geo_Lshaped_two_patches.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_diff  = @(x, y) ones(size(x));
%! k = 1; % Constant that characterizes the singularity
%! problem_data.f = @(x, y) zeros (size (x));
%! problem_data.h = @(x, y, ind) singular_function_laplace (x, y, k);
%! problem_data.uex     = @(x, y) singular_function_laplace (x, y, k);
%! problem_data.graduex = @(x, y) singular_function_maxwell (x, y, k);
%! method_data.degree      = [3 3];        % Degree of the splines
%! method_data.regularity  = [2 2];        % Regularity of the splines
%! method_data.nsub_coarse = [2 2];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
%! method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
%! method_data.nquad       = [4 4];        % Points for the Gaussian quadrature rule
%! method_data.space_type  = 'simplified'; % 'simplified' (only children functions) or 'standard' (full basis)
%! method_data.truncated   = 0;            % 0: False, 1: True
%! adaptivity_data.flag = 'functions';
%! adaptivity_data.C0_est = 1.0;
%! adaptivity_data.mark_param = .5;
%! adaptivity_data.mark_strategy = 'MS';
%! adaptivity_data.max_level = 10;
%! adaptivity_data.max_ndof = 15000;
%! adaptivity_data.num_max_iter = 5;
%! adaptivity_data.max_nel = 15000;
%! adaptivity_data.tol = 1e-10;
%! plot_data.print_info = false;
%! plot_data.plot_hmesh = false;
%! plot_data.plot_discrete_sol = false;
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.iter, 5)
%! assert (solution_data.ndof, [45 61 83 104 149])
%! assert (solution_data.nel, [8 20 38 62 98])
%! assert (solution_data.err_h1s, [1.089825381386800e-01 7.377184523296307e-02 4.889161568362813e-02 3.25822537253584e-02 2.133745968166025e-02], 1e-15)
