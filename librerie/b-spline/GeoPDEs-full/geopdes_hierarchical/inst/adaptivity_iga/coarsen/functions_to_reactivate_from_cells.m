% FUNCTIONS_TO_REACTIVATE_FROM_CELLS: compute the indices of the active
%  functions that have to be activated during coarsening.
%
%   fun_indices = functions_to_reactivate_from_cells (hspace, hmsh, marked)
%
% INPUT:
%
%   hspace: object representing the (coarse) space of hierarchical splines (see hierarchical_space)
%   hmsh:   object representing the refined hierarchical mesh (see hierarchical_mesh)
%   marked: cell array with the indices, in the Cartesian grid, of the elements that have to 
%             be reactivated for each level
%
% OUTPUT:
%   fun_indices: cell array with the indices of functions to be removed, for
%       each level, in the numbering of the tensor product space
%
% Copyright (C) 2016 Eduardo M. Garau, Rafael Vazquez
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

function fun_indices = functions_to_reactivate_from_cells (hmsh, hspace, M)

fun_indices = cell (hspace.nlevels,1);

for lev = 1:hspace.nlevels-1
  if (~isempty(M{lev}))
    fun_indices{lev} = union (fun_indices{lev}, sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), M{lev}));
    fun_indices{lev} = intersect (hspace.deactivated{lev}, fun_indices{lev});
    
    if (strcmpi (hspace.type, 'simplified'))
      children = intersect (hspace_get_children (hspace, lev, fun_indices{lev}), hspace.deactivated{lev+1});
      aux_deact = setdiff (hspace.deactivated{lev}, fun_indices{lev});
      for ifun = children(:)'
        if (~isempty (intersect (hspace_get_parents (hspace, lev+1, ifun), aux_deact)))
          children = setdiff (children, ifun);
        end
      end
      fun_indices{lev+1} = children;
    end
  end
end

end
