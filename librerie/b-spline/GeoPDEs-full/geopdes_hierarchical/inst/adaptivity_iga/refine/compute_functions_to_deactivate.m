% COMPUTE_FUNCTIONS_TO_DEACTIVATE: compute the indices of the active functions that have to be deactivated.
%
%   fun_indices = compute_functions_to_deactivate (hspace, hmsh, marked, flag)
%
% INPUT:
%
%   hspace: object representing the (coarse) space of hierarchical splines (see hierarchical_space)
%   hmsh:   object representing the refined hierarchical mesh (see hierarchical_mesh)
%   marked: cell array with the indices, in the tensor product space, of the marked 
%             element/functions for each level
%   flag:   the refinement strategy, marking either 'elements' or 'functions'
%
% OUTPUT:
%   fun_indices: cell array with the indices of functions to be deactivated, for
%       each level, in the numbering of the tensor product space
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

function fun_indices = compute_functions_to_deactivate (hmsh, hspace, M, flag)

fun_indices = cell (hspace.nlevels,1);

for lev = 1:hspace.nlevels
  if (~isempty(M{lev}))
  % Computation of candidate functions to be deactivated. Only active
  %  functions which are not already marked are considered
    switch (flag)
      case ('functions')
        fun_indices{lev} = sp_get_neighbors (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), M{lev});
        fun_indices{lev} = setdiff (intersect (fun_indices{lev}, hspace.active{lev}), M{lev});
      case ('elements')
        fun_indices{lev} = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), M{lev});
        fun_indices{lev} = intersect (fun_indices{lev}, hspace.active{lev});
    end

  % Computation of functions that in fact have to be deactivated
    [~, cells_per_fun] = sp_get_cells (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), fun_indices{lev});
    flag_ell = cellfun (@(x) isempty (intersect (x, hmsh.active{lev})), cells_per_fun);
    fun_indices{lev} = fun_indices{lev}(flag_ell == 1);
        
    if (strcmpi (flag,'functions'))
      fun_indices{lev} = union (fun_indices{lev}, M{lev});
    end
    fun_indices{lev} = fun_indices{lev}(:);
  end
end

end