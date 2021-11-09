% PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_cube.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));
problem_data.grad_c_diff = @(x, y, z) cat (1, ...
            reshape (zeros(size(x)), [1, size(x)]), ...
            reshape (zeros(size(x)), [1, size(x)]), ...
            reshape (zeros(size(x)), [1, size(x)]));

% Source and boundary terms
C = 100;
normax2 = @(x,y,z) ((x-.5).^2+(y-.5).^2+(z-.5).^2);
uex = @(x,y,z) exp(-C*normax2(x,y,z));
problem_data.f = @(x, y, z) 4*C*(6/4-C*normax2(x,y,z)).*uex(x,y,z);
problem_data.h = @(x, y, z, ind) uex(x,y,z);

% Exact solution (optional)
problem_data.uex     = uex;
problem_data.graduex = @(x, y, z) -2*C*cat (1, ...
            reshape (uex(x,y,z).*(x-.5), [1, size(x)]), ...
            reshape (uex(x,y,z).*(y-.5), [1, size(x)]), ...
            reshape (uex(x,y,z).*(z-.5), [1, size(x)]));

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [2 2 2];      % Degree of the splines
method_data.regularity  = [1 1 1];      % Regularity of the splines
method_data.nsub_coarse = [2 2 2];      % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2 2];      % Number of subdivisions for each refinement
method_data.nquad       = [3 3 3];      % Points for the Gaussian quadrature rule
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
adaptivity_data.max_ndof = 5000;
adaptivity_data.num_max_iter = 5;
adaptivity_data.max_nel = 5000;
adaptivity_data.tol = 1e-5;

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

[geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace(problem_data, method_data, adaptivity_data, plot_data);

% EXPORT VTK FILE
npts = [51 51 51];
output_file = 'laplace_adaptivity_cube.vts';
sp_to_vtk (u, hspace, geometry, npts, output_file, {'solution', 'gradient', 'laplacian'}, {'value', 'gradient', 'laplacian'})


%!test
%! problem_data.geo_name = 'geo_cube.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4 5 6];
%! problem_data.c_diff  = @(x, y, z) ones(size(x));
%! C = 100;
%! normax2 = @(x,y,z) ((x-.5).^2+(y-.5).^2+(z-.5).^2);
%! uex = @(x,y,z) exp(-C*normax2(x,y,z));
%! problem_data.f = @(x, y, z) 4*C*(6/4-C*normax2(x,y,z)).*uex(x,y,z);
%! problem_data.h = @(x, y, z, ind) uex(x,y,z);
%! problem_data.uex     = uex;
%! problem_data.graduex = @(x, y, z) -2*C*cat (1, ...
%!             reshape (uex(x,y,z).*(x-.5), [1, size(x)]), ...
%!             reshape (uex(x,y,z).*(y-.5), [1, size(x)]), ...
%!             reshape (uex(x,y,z).*(z-.5), [1, size(x)]));
%! method_data.degree      = [2 2 2];      % Degree of the splines
%! method_data.regularity  = [1 1 1];      % Regularity of the splines
%! method_data.nsub_coarse = [2 2 2];      % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
%! method_data.nsub_refine = [2 2 2];      % Number of subdivisions for each refinement
%! method_data.nquad       = [3 3 3];      % Points for the Gaussian quadrature rule
%! method_data.space_type  = 'simplified'; % 'simplified' (only children functions) or 'standard' (full basis)
%! method_data.truncated   = 0;            % 0: False, 1: True
%! adaptivity_data.flag = 'functions';
%! adaptivity_data.C0_est = 1.0;
%! adaptivity_data.mark_param = .5;
%! adaptivity_data.mark_strategy = 'MS';
%! adaptivity_data.max_level = 10;
%! adaptivity_data.max_ndof = 5000;
%! adaptivity_data.num_max_iter = 4;
%! adaptivity_data.max_nel = 5000;
%! adaptivity_data.tol = 1e-5;
%! plot_data.print_info = false;
%! plot_data.plot_hmesh = false;
%! plot_data.plot_discrete_sol = false;
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace(problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.iter, 4)
%! assert (solution_data.ndof, [64 216 1000 1208])
%! assert (solution_data.nel, [8 64 512 960])
%! assert (solution_data.err_h1s, [1.291487658740688e+00 6.708611152128069e-01 4.533759941100940e-01 6.012296757586015e-02], 1e-15)
