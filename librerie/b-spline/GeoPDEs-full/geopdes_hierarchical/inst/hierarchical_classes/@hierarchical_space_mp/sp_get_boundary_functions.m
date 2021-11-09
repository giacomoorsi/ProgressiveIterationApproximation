% SP_GET_BOUNDARY_FUNCTIONS: indices of the degrees of freedom of the given boundaries.
%
%   dofs = sp_get_boundary_functions (hspace, hmsh, sides)
%
% INPUT:
%
%  hspace: object representing the hierarchical space of trial functions (see hierarchical_space_mp)
%  hmsh:   object representing the hierarchical mesh (see hierarchical_mesh_mp)
%  sides:  boundary sides from which we want to compute the indices
%
% OUTPUT:
%
%  dofs: global numbering of the boundary basis functions
%
% Copyright (C) 2019 Rafael Vazquez
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

function dofs = sp_get_boundary_functions (hspace, hmsh, drchlt_sides)

  boundaries = hmsh.mesh_of_level(1).boundaries;
  Nbnd = cumsum ([0, boundaries.nsides]);
  bnd_dofs = [];

  Nf = cumsum ([0, hspace.boundary.ndof_per_level]);
  for iref = drchlt_sides
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    for lev = 1:hspace.boundary.nlevels
%       boundary_gnum = hspace.boundary.space_of_level(lev).gnum;
      boundary_gnum = cellfun (@(x) x(:), hspace.boundary.space_of_level(lev).gnum, 'UniformOutput', false);
      global_boundary_dofs = vertcat (boundary_gnum{iref_patch_list});
      [~,bnd_dofs_lev] = intersect (hspace.boundary.active{lev}, global_boundary_dofs(:));
      bnd_dofs = union (bnd_dofs, Nf(lev) + bnd_dofs_lev);
    end
  end
  
  dofs = hspace.boundary.dofs(bnd_dofs);

end
