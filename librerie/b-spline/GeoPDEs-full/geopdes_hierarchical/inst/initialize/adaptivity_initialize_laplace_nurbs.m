% ADAPTIVITY_INITIALIZE_LAPLACE: initialize a hierarchical mesh and a hierarchical space with one-single level.
%
% [hmsh, hspace] = adaptivity_initialize_laplace (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. For this function, it must contain the fields:
%    - geo_name:     name of the file containing the geometry
%
%  method_data : a structure with discretization data. For this function, it contains the fields:
%    - degree:      degree of the spline functions.
%    - regularity:  continuity of the spline functions.
%    - nsub_coarse: number of subelements with respect to the geometry mesh (1 leaves the mesh unchanged)
%    - nsub_refine: number of subelements to be added at each refinement step (2 for dyadic)
%    - nquad:       number of points for Gaussian quadrature rule
%    - space_type:  'simplified' (only children of removed functions) or 'standard' (full hierarchical basis)
%    - truncated:   false (classical basis) or true (truncated basis)
%
% OUTPUT:
%    hmsh:     object representing the hierarchical mesh (see hierarchical_mesh)
%    hspace:   object representing the space of hierarchical splines (see hierarchical_space)
%    geometry: geometry structure (see geo_load)
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function [hmsh, hspace, geometry] = adaptivity_initialize_laplace_nurbs (problem_data, method_data)

[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (problem_data.geo_name);
geometry = rmfield (geometry, 'name');

npatch = numel (geometry);
for iptc = 1:npatch
  degelev  = max (method_data.degree - (geometry(iptc).nurbs.order-1), 0);
  nurbs    = nrbdegelev (geometry(iptc).nurbs, degelev);
  [~, zeta, nknots] = kntrefine (nurbs.knots, method_data.nsub_coarse-1, nurbs.order-1, method_data.regularity);

  nurbs = nrbkntins (nurbs, nknots);
  geometry(iptc) = orderfields (geo_load (nurbs), geometry);
    
  rule     = msh_gauss_nodes (method_data.nquad);
  [qn, qw] = msh_set_quad_nodes (zeta, rule);
  msh{iptc}   = msh_cartesian (zeta, qn, qw, geometry(iptc));
  space{iptc} = sp_nurbs (nurbs, msh{iptc});
end

if (npatch == 1)
  msh = msh{1};
  space = space{1};
  hmsh     = hierarchical_mesh (msh, method_data.nsub_refine);
  hspace   = hierarchical_space (hmsh, space, method_data.space_type, method_data.truncated);
else
  msh = msh_multipatch (msh, boundaries);
  space = sp_multipatch (space, msh, interfaces, boundary_interfaces);
  hmsh     = hierarchical_mesh_mp (msh, method_data.nsub_refine);
  hspace   = hierarchical_space_mp (hmsh, space, method_data.space_type, method_data.truncated);
end


end