% SP_DRCHLT_L2_PROJ: assign the degrees of freedom of Dirichlet boundaries through an L2 projection.
%
%   [u, dofs] = sp_drchlt_l2_proj (hspace, hmsh, h, sides)
%
% INPUT:
%
%  hspace: object representing the hierarchical space of trial functions (see hierarchical_space_mp)
%  hmsh:   object representing the hierarchical mesh (see hierarchical_mesh_mp)
%  h:      function handle to compute the Dirichlet condition
%  sides:  boundary sides on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: global numbering of the corresponding basis functions
%
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
% Copyright (C) 2015 Eduardo M. Garau
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

function [u, dofs] = sp_drchlt_l2_proj (hspace, hmsh, h, drchlt_sides)

  M = spalloc (hspace.boundary.ndof, hspace.boundary.ndof, 3*hspace.boundary.ndof);
  rhs  = zeros (hspace.boundary.ndof, 1);

  boundaries = hmsh.mesh_of_level(1).boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  bnd_dofs = [];

  Nf = cumsum ([0, hspace.boundary.ndof_per_level]);
  for iref = drchlt_sides
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    href = @(varargin) h(varargin{:}, iref);
    f_one = @(varargin) ones (size(varargin{1}));
    
    M = M + op_u_v_hier (hspace.boundary, hspace.boundary, hmsh.boundary, f_one, iref_patch_list);
    rhs = rhs + op_f_v_hier (hspace.boundary, hmsh.boundary, href, iref_patch_list);

    for lev = 1:hspace.boundary.nlevels
%       boundary_gnum = hspace.boundary.space_of_level(lev).gnum;
      boundary_gnum = cellfun (@(x) x(:), hspace.boundary.space_of_level(lev).gnum, 'UniformOutput', false);
      global_boundary_dofs = vertcat (boundary_gnum{iref_patch_list});
      [~,bnd_dofs_lev] = intersect (hspace.boundary.active{lev}, global_boundary_dofs(:));
      bnd_dofs = union (bnd_dofs, Nf(lev) + bnd_dofs_lev);
    end
  end
  
  u = M(bnd_dofs,bnd_dofs) \ rhs(bnd_dofs, 1);
  dofs = hspace.boundary.dofs(bnd_dofs);

end
