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
% Choose 1.5 < Cx, Cy < 2.5. Thus, the solution belongs H2\H3
Cx = 1.6;
Cy = 2.4;
problem_data.f = @(x,y) Cx*((Cx+1)*x.^(Cx-1)-(Cx-1)*x.^(Cx-2)).*y.^Cy.*(1-y)+...
            Cy*((Cy+1)*y.^(Cy-1)-(Cy-1)*y.^(Cy-2)).*x.^Cx.*(1-x);
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));

% Exact solution (optional)
problem_data.uex = @(x,y) x.^Cx.*(1-x).*y.^Cy.*(1-y);
problem_data.graduex = @(x,y) cat (1, ...
            reshape ((Cx*x.^(Cx-1).*(1-x)-x.^Cx).*y.^Cy.*(1-y), [1, size(x)]), ...
            reshape ((Cy*y.^(Cy-1).*(1-y)-y.^Cy).*x.^Cx.*(1-x), [1, size(x)]));
            
        
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
adaptivity_data.max_level = 15;
adaptivity_data.max_ndof = 10000;
adaptivity_data.num_max_iter = 10;
adaptivity_data.max_nel = 10000;
adaptivity_data.tol = 1e-5;

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

[geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);

% EXPORT VTK FILE
npts = [51 51];
output_file = 'laplace_adaptivity_square_ex3.vts';
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
%! Cx = 1.6;
%! Cy = 2.4;
%! problem_data.f = @(x,y) Cx*((Cx+1)*x.^(Cx-1)-(Cx-1)*x.^(Cx-2)).*y.^Cy.*(1-y)+...
%!             Cy*((Cy+1)*y.^(Cy-1)-(Cy-1)*y.^(Cy-2)).*x.^Cx.*(1-x);
%! problem_data.g = @(x, y, ind) zeros(size(x));
%! problem_data.h = @(x, y, ind) zeros(size(x));
%! problem_data.uex = @(x,y) x.^Cx.*(1-x).*y.^Cy.*(1-y);
%! problem_data.graduex = @(x,y) cat (1, ...
%!             reshape ((Cx*x.^(Cx-1).*(1-x)-x.^Cx).*y.^Cy.*(1-y), [1, size(x)]), ...
%!             reshape ((Cy*y.^(Cy-1).*(1-y)-y.^Cy).*x.^Cx.*(1-x), [1, size(x)]));
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
%! adaptivity_data.max_level = 15;
%! adaptivity_data.max_ndof = 10000;
%! adaptivity_data.num_max_iter = 5;
%! adaptivity_data.max_nel = 10000;
%! adaptivity_data.tol = 1e-5;
%! plot_data.print_info = false;
%! plot_data.plot_hmesh = false;
%! plot_data.plot_discrete_sol = false;
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.iter, 5)
%! assert (solution_data.ndof, [25 49 73 109 178]);
%! assert (solution_data.nel, [4 16 37 67 130]);
%! assert (solution_data.err_h1s, [1.585089393282615e-03 5.245888575707342e-04 2.975940577220350e-04 1.490977135010352e-04 9.225171525608927e-05], 1e-15)
%! assert (solution_data.flag, 2)
%!
%! method_data.space_type  = 'standard';
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.iter, 5)
%! assert (solution_data.ndof, [25 49 73 109 178]);
%! assert (solution_data.nel, [4 16 37 67 130]);
%! assert (solution_data.err_h1s, [1.585089393282615e-03 5.245888575707342e-04 2.975940577220350e-04 1.490977135010352e-04 9.225171525608927e-05], 1e-15)
%! assert (solution_data.flag, 2)
%!
%! adaptivity_data.flag = 'elements';
%! [geometry, hmsh, hspace, u, solution_data] = adaptivity_laplace (problem_data, method_data, adaptivity_data, plot_data);
%! assert (solution_data.iter, 5)
%! assert (solution_data.ndof, [25 40 58 79 109]);
%! assert (solution_data.nel, [4 13 25 49 73]);
%! assert (solution_data.err_h1s, [1.585089393282615e-03 6.641906171007896e-04 3.337431928336381e-04 2.239401521718252e-04 1.475543816385030e-04], 1e-15)
%! assert (solution_data.flag, 2)
