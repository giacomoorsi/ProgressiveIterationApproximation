% HIERARCHICAL_SPACE: constructor of the class for hierarchical spaces.
%
%    function hspace = hierarchical_space (hmsh, space, [space_type, truncated, regularity])
%
% INPUT
%    hmsh:       an object of the class hierarchical_mesh (see hierarchical_mesh)
%    space:      the coarsest space, an object of the class sp_scalar (see sp_scalar or sp_vector)
%    space_type: select which kind of hierarchical space to construct. The options are
%                - 'standard',   the usual hierachical splines space (default value)
%                - 'simplified', a simplified basis, were only children of removed functions are activated
%    truncated:  decide whether the basis will be truncated or not (not truncated by default)
%    regularity: will be used for refinement. For vectors, it should be given in a cell array. By default it is degree minus one
%
% OUTPUT:
%    hspace: hierarchical_space object, which contains the following fields and methods
%     (fields between brackets are only used for vector-valued spaces)
% 
%    FIELD_NAME     TYPE                    DESCRIPTION
%    ncomp          (scalar)                number of components of the space
%    [ncomp_param]  (scalar)                number of components of the functions of the space in the parametric domain (usually equal to msh.ndim)
%    type           (string)                'standard' or 'simplified'
%    ndof           (scalar)                total number of active functions 
%    [comp_dofs]    (1 x ncomp_param cell array) indices of the degrees of freedom for each component
%    nlevels        (scalar)                the number of levels
%    space_of_level (1 x nlevels)           tensor product space of each level, with 1d evaluations on the mesh of the same level (see sp_bspline)
%    Proj           (hmsh.nlevels-1 x 1 cell-array) each level contains a second cell array
%      (for level)  (ncomp_param x ndim cell-array) 
%                                           the coefficients relating 1D splines of two consecutive levels
%                                           Proj{l}{i} is a matrix of size N_{l+1} x N_l where N_l is the number 
%                                           of univariate functions of level l in the direction l, such that
%                                           a function B_{k,l} = \sum_j c^k_j B_{j,l+1}, and c^k_j = Proj{l}{i}(j,k)
%                                           For vectors, it takes the form Proj{l}{c,i}, 
%                                           where c is the component in the parametric domain
%    ndof_per_level (1 x nlevels array)     number of active functions on each level
%    active        (1 x nlevels cell-array) List of active functions on each level
%    coeff_pou     (ndof x 1)               coefficientes to form the partition of the unity in the hierarchical space
%    deactivated   (1 x nlevels cell-array) List of deactivated functions on each level
%    Csub          (1 x hmsh.nlevels cell-array) Sparse matrices for changing basis. For each level, represent active functions of previous levels
%                                            as linear combinations of splines (active and inactive) of the current level
%    Csub_row_indices (1 x hmsh.nlevels cell-array) indices of the rows stored in Csub. 
%                                            This allows to save memory space.
%    boundary      (2 x ndim array)         a hierarchical space representing the restriction to the boundary
%    dofs          (1 x ndof array)         only for boundary spaces, degrees of freedom that do not vanish on the boundary
%    regularity    (1 x ndim array)         the regularity of the space, used during refinement to add a new level 
%               or (1 x ncomp cell-array)
%
%    METHODS
%    Methods for post-processing, which require a computed vector of degrees of freedom
%      sp_to_vtk:             export the solution to a VTK file, in a structured grid of points
%      sp_eval:               evaluate the solution in a Cartesian grid of points
%      sp_plot_solution:      plot the computed solution, given the degrees of freedom
%      hspace_eval_hmsh:      evaluate the solution in the quadrature points of the corresponding hierarchical mesh
%      sp_l2_error:           compute the error in L2 norm
%      sp_h1_error:           compute the error in H1 norm
%
%    Methods for matrix and vector assembly, and boundary conditions
%      op_gradu_gradv_hier:   assemble the stiffness matrix
%      op_u_v_hier:           assemble the mass matrix
%      op_f_v_hier:           assemble the right-hand side
%      sp_drchlt_l2_proj:     compute the boundary degrees of freedom using the L2-projection
%
%    Methods to obtain a relation between different levels
%      hspace_get_children:   get the children of a given set of functions
%      hspace_get_parents:    get the parents of a given set of functions
%      hspace_subdivision_matrix: compute the matrix of basis change (the same stored in Csub)
%
%    Other methods
%      hspace_refine:         refine the hierarchical space
%      hspace_coarsen:        coarsen the hierarchical space
%      hspace_add_new_level:  add a new level, initialized without active functions
%      hspace_remove_empty_levels: remove the finest level, if it is empty
%      hspace_in_finer_mesh:  compute the same space in a finer hierarchical mesh
%      hspace_admissibility_class: check the admissibility class of the associated mesh.
%      hspace_check_partition_of_unity: check whether the computed coefficients
%                             for the partition of unity are correct (used for debugging)
%      sp_get_boundary_functions: get the degrees of freedom of a given boundary
%
% For details about the 'simplified' hierarchical space:
%    A. Buffa, E. M. Garau, Refinable spaces and local approximation estimates 
%     for hierarchical splines, IMA J. Numer. Anal., (2016)
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2017-2019 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function hspace = hierarchical_space (hmsh, space, varargin)

if (isa (space, 'sp_scalar'))
  is_scalar = true;
  regularity = space.degree - 1;
elseif (isa (space, 'sp_vector'))
  is_scalar = false;
  for icomp = 1:space.ncomp_param
    regularity{icomp} = space.scalar_spaces{icomp}.degree-1;
  end
else
  error ('Unknown space type')
end

default_values = {'standard', false, regularity};
default_values(1:numel(varargin)) = varargin;
[space_type, truncated, regularity] = default_values{:};

hspace.ncomp = space.ncomp;
hspace.type = space_type;
hspace.truncated = truncated;

hspace.nlevels = 1;
hspace.ndof = space.ndof;
hspace.ndof_per_level = space.ndof;
hspace.space_of_level = space;

hspace.active{1} = (1:space.ndof)';
hspace.deactivated{1} = [];

hspace.coeff_pou = ones (space.ndof, 1);
if (is_scalar)
  hspace.Proj = cell (0, 1);
  hspace.ncomp_param = 1;
  hspace.comp_dofs = [];
else
  hspace.Proj = cell (0, 1);
  hspace.ncomp_param = space.ncomp_param;
  aux = 0;
  for icomp = 1:space.ncomp_param
    hspace.comp_dofs{icomp} = aux+(1:space.scalar_spaces{icomp}.ndof);
    aux = aux + space.scalar_spaces{icomp}.ndof;
  end
end
hspace.Csub{1} = speye (space.ndof);
hspace.Csub_row_indices{1} = 1:space.ndof;

hspace.dofs = [];
hspace.adjacent_dofs = [];

transform = hspace.space_of_level(1).transform;
if (~isempty (hmsh.boundary) && ~isempty (space.boundary))
  if (hmsh.ndim > 1)
    for iside = 1:numel (hmsh.boundary)
% This check is necessary for the regularity
      ind = setdiff (1:hmsh.ndim, ceil(iside/2));
      if (is_scalar)
        bnd_regularity = regularity(ind);
      else
        if (strcmpi (transform, 'grad-preserving'))
          comps = 1:space.ncomp;
          bnd_regularity = cellfun (@(x) x(ind), regularity(comps), 'UniformOutput', false);
        elseif (strcmpi (transform, 'curl-preserving'))
          comps = setdiff (1:hmsh.ndim, ceil(iside/2)); % comps =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %comps = [2 2 1 1] in 2D;
          bnd_regularity = cellfun (@(x) x(ind), regularity(comps), 'UniformOutput', false);
%         elseif (strcmpi (transform, 'div-preserving'))
%           comps = ceil (iside/2); % comps =[1, 1, 2, 2, 3, 3] in 3D, %comps = [1, 1, 2, 2] in 2D;
%           bnd_regularity = regularity{comps}(ind);
        end
      end

      boundary = hierarchical_space (hmsh.boundary(iside), space.boundary(iside), space_type, truncated, bnd_regularity);
      boundary.dofs = space.boundary(iside).dofs;
      if (is_scalar)
        boundary.adjacent_dofs = space.boundary(iside).adjacent_dofs;
      else
        boundary.adjacent_dofs = [];
      end
      hspace.boundary(iside) = boundary;
    end
  elseif (hmsh.ndim == 1)
    hspace.boundary(1).dofs = space.boundary(1).dofs;
    hspace.boundary(2).dofs = space.boundary(2).dofs;
  end
else
  hspace.boundary = [];
end

hspace.regularity = regularity;

hspace = class (hspace, 'hierarchical_space');

end
