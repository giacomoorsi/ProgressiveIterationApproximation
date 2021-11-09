% HMSH_EVAL_BOUNDARY_SIDE: evaluate the parameterization in one boundary side of the domain.
%
%     msh_side = hmsh_eval_boundary_side (hmsh, iside, [element_list]);
%
% INPUTS:
%     
%    msh:   mesh object (see hierarchical_mesh)
%    iside: number of the boundary side to compute, from 1 to 2*msh.ndim (see the file geo_specs for documentation about face numbering)
%    element_list: elements on which to compute the parametrization. They
%      may be from different levels. All elements are taken, by default.
%
% OUTPUT:
%
%     msh_side: structure that contains the following fields
%
%     FIELD_NAME    (SIZE)                  DESCRIPTION
%     side_number   (scalar)                  number of the side
%     nel           (scalar)                  number of elements of the boundary side
%     nel_dir       (1 x ndim vector)         number of elements in each parametric direction
%     nqn           (scalar)                  number of quadrature nodes per element
%     nqn_dir       (1 x ndim vector)         number of quadrature nodes per element in each parametric direction
%     quad_nodes    (ndim x nqn x nel vector) coordinates of the quadrature nodes in parametric domain
%     quad_weights  (nqn x nel vector)        weights associated to the quadrature nodes
%     geo_map       (rdim x nqn x nel vector) physical coordinates of the quadrature nodes
%     geo_map_jac   (rdim x ndim x nqn x nel) Jacobian matrix of the map evaluated at the quadrature nodes
%     jacdet        (nqn x nel)               element of length, area, volume (if rdim = ndim, determinant of the Jacobian)
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2014, 2015, 2017 Rafael Vazquez
% Copyright (C) 2017 Luca Coradello
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

function msh_side = hmsh_eval_boundary_side (hmsh, iside, element_list)

hmsh_bnd = hmsh.boundary(iside);
if (nargin == 2)
  element_list = 1:hmsh_bnd.nel;
end

if (any (element_list > hmsh_bnd.nel))
  error ('Element index greater than the number of elements')
end

msh_side.side_number = iside;
msh_side.ndim = hmsh_bnd.ndim;
msh_side.rdim = hmsh_bnd.rdim;
msh_side.nel = numel (element_list);

nqn = hmsh_bnd.mesh_of_level(1).nqn;
if (isempty (element_list))
  msh_side.quad_weights = [];
  msh_side.geo_map = [];
  msh_side.geo_map_jac = [];
  msh_side.geo_map_der2 = [];
  msh_side.jacdet = [];
  msh_side.element_size = [];
  return
else
  nel = numel (element_list);
  msh_side.nqn = nqn;
  msh_side.quad_weights = zeros (nqn, nel);
  msh_side.geo_map = zeros (hmsh_bnd.rdim, nqn, nel);
  msh_side.geo_map_jac = zeros (hmsh_bnd.rdim, hmsh_bnd.ndim, nqn, nel);
  msh_side.geo_map_der2 = zeros (hmsh_bnd.rdim, hmsh_bnd.ndim, hmsh_bnd.ndim, nqn, nel);
  msh_side.jacdet = zeros (nqn, nel);
  msh_side.element_size = zeros (1, nel);
  msh_side.normal = zeros (hmsh_bnd.rdim, nqn, nel);
end


% This avoids an error with subsasgn. Only the side "iside" will be used.
boundary = hmsh_bnd.mesh_of_level(1);
for side = 2:iside
  boundary(side) = hmsh_bnd.mesh_of_level(1);
end

shifting_indices = cumsum ([0 hmsh_bnd.nel_per_level]);
for ilev = 1:hmsh_bnd.nlevels
  elements_of_level = shifting_indices(ilev)+1:shifting_indices(ilev+1);
  [~,level_indices,input_indices] = intersect (elements_of_level, element_list);
  if (~isempty (level_indices))
    mesh_of_level = hmsh.mesh_of_level(ilev);
    boundary(iside) = hmsh_bnd.mesh_of_level(ilev);
    mesh_of_level.boundary = boundary;
    msh_aux = msh_eval_boundary_side (mesh_of_level, iside, hmsh_bnd.active{ilev}(level_indices));

    msh_side.quad_weights(:,input_indices) = msh_aux.quad_weights;
    msh_side.geo_map(:,:,input_indices) = msh_aux.geo_map;
    msh_side.geo_map_jac(:,:,:,input_indices) = msh_aux.geo_map_jac;
    msh_side.geo_map_der2(:,:,:,:,input_indices) = msh_aux.geo_map_der2;
    msh_side.jacdet(:,input_indices) = msh_aux.jacdet;
    msh_side.element_size(:,input_indices) = msh_aux.element_size;
    msh_side.normal(:,:,input_indices) = msh_aux.normal;
  end
end

end
