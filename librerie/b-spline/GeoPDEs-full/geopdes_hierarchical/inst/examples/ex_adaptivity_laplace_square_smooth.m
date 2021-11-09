% PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));
problem_data.grad_c_diff = @(x, y) cat (1, ...
            reshape (zeros(size(x)), [1, size(x)]), ...
            reshape (zeros(size(x)), [1, size(x)]));
        
% Source and boundary terms
problem_data.f = @(x, y) zeros (size (x));
problem_data.g = @test_square_g_nmnn;
problem_data.h = @(x, y, ind) exp (x) .* sin(y);

% Exact solution (optional)
problem_data.uex     = @(x, y) exp (x) .* sin (y);
problem_data.graduex = @(x, y) cat (1, ...
                       reshape (exp(x).*sin(y), [1, size(x)]), ...
                       reshape (exp(x).*cos(y), [1, size(x)]));

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = [3 3];        % Degree of the splines
method_data.regularity  = [2 2];        % Regularity of the splines
method_data.nsub_coarse = [3 3];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
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
adaptivity_data.num_max_iter = 7;
adaptivity_data.max_nel = 15000;
adaptivity_data.tol = 1e-10;

% GRAPHICS
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;
plot_data.print_info = true;

[geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);

% EXPORT VTK FILE
npts = [51 51];
output_file = 'laplace_adaptivity_square_ex1.vts';
sp_to_vtk (u, hspace, geometry, npts, output_file, {'solution', 'gradient', 'laplacian'}, {'value', 'gradient', 'laplacian'})

% Plot in Octave/Matlab
[eu, F] = sp_eval (u, hspace, geometry, npts);
figure; subplot (1,2,1)
surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu)
subplot(1,2,2)
surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze (problem_data.uex(F(1,:,:), F(2,:,:))));

%!test
%! problem_data.geo_name = 'geo_square.txt';
%! problem_data.nmnn_sides   = [];
%! problem_data.drchlt_sides = [1 2 3 4];
%! problem_data.c_diff  = @(x, y) ones(size(x));
%! problem_data.grad_c_diff = @(x, y) zeros ([2, size(x)]);
%! problem_data.f = @(x, y) zeros (size (x));
%! problem_data.g = @test_square_g_nmnn;
%! problem_data.h = @(x, y, ind) exp (x) .* sin(y);
%! problem_data.uex     = @(x, y) exp (x) .* sin (y);
%! problem_data.graduex = @(x, y) cat (1, ...
%!                        reshape (exp(x).*sin(y), [1, size(x)]), ...
%!                        reshape (exp(x).*cos(y), [1, size(x)]));
%! method_data.degree      = [3 3];        % Degree of the splines
%! method_data.regularity  = [2 2];        % Regularity of the splines
%! method_data.nsub_coarse = [3 3];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
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
%! assert (solution_data.ndof, [36 81 210 372 615]);
%! assert (solution_data.nel, [9 36 138 291 510]);
%! assert (solution_data.err_h1s, [2.448941825078764e-04, 3.245165315459014e-05, 4.555200816633961e-06, ...
%!           2.449082638340557e-06, 6.648399812699606e-07], 1e-15)
%! assert (solution_data.flag, 2)
%!
%! adaptivity_data.flag = 'elements';
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.iter, 5)
%! assert (solution_data.ndof, [36 45 60 72 84]);
%! assert (solution_data.nel, [9 18 27 33 42]);
%! assert (solution_data.err_h1s, [2.448941825078764e-04, 1.959352092852789e-04, 8.135128968531835e-05, 4.519709525363875e-05, 3.231472797704601e-05], 1e-15)
%! assert (solution_data.flag, 2)
%!
%! method_data.space_type  = 'standard';
%! method_data.truncated = true;
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.iter, 5)
%! assert (solution_data.ndof, [36 45 60 72 84]);
%! assert (solution_data.nel, [9 18 27 33 42]);
%! assert (solution_data.err_h1s, [2.448941825078764e-04, 1.959352092852789e-04, 8.135128968531835e-05, 4.519709525363875e-05, 3.231472797704601e-05], 1e-15)
%! assert (solution_data.flag, 2)
