% COMPUTE_CELLS_TO_REFINE: compute the indices of the cells that have to be refined, when 
%  the refinement strategy is based on marking functions.
%
%   [marked_elem, indices] = compute_cells_to_refine (hspace, hmsh, marked)
%
% INPUT:
%
%   hspace: object representing the space of hierarchical splines (see hierarchical_space)
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   marked_fun: cell array with the indices, in the tensor product space, of the marked functions
%            for each level
%
% OUTPUT:
%
%   marked_elem: cell array with the indices, in the tensor Cartesian grid, of the elements
%                 to be refined for each level
%   indices:     relative position of the marked elements in the numbering of the hierarchical mesh
%                 within the level, that is, in hmsh.active{lev}
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

function [ME, ind] = compute_cells_to_refine (hspace, hmsh, MF)

ME = cell (hmsh.nlevels, 1);
ind = cell (hmsh.nlevels, 1);

for lev = 1:hspace.nlevels
  if (~isempty(MF{lev}))
    ME{lev} = sp_get_cells (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), MF{lev});
    [ME{lev}, ~, ind{lev}] = intersect (ME{lev}, hmsh.active{lev});
  end
end

end