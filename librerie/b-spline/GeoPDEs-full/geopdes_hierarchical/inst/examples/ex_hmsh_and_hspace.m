% ex_hmsh_and_hspace
%
% Test file. It uses build_hspace_from_cells to generate a hierarchical mesh
%  and space from the set of refined cells of each level

dim = 2;    % number of parametric directions
n = 6;      % number of initial elements in each direction, [6 6] is also valid
degree = 2; % polynomial degree, [2 2] is also valid.
simplified = false; % Use the reduced (simplified) space of HB-splines.
truncated = true;   % Use the truncated space or not.

% Give either the list of 'refined' cells, or the list of 'active' cells
cell_type = 'refined';
cells{1} = [15 16 21 22];
cells{2} = [53 54];

[hmsh,hspace] = build_hspace_from_cells (dim, degree, n, cells, cell_type, simplified, truncated);
hmsh_plot_cells (hmsh)