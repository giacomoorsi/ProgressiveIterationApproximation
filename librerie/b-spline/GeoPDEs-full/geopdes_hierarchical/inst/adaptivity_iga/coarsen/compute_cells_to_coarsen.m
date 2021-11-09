% COMPUTE_CELLS_TO_COARSEN: compute the indices of the cells that have to be coarsened, when 
%  the refinement strategy is based on marking functions.
%
%   [marked_elem, indices] = compute_cells_to_coarsen (hspace, hmsh, marked)
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
%                 to be coarsened for each level
%   indices:     relative position of the marked elements in the numbering of the hierarchical mesh
%                 within the level, that is, in hmsh.active{lev}
%
% From the list of (active) marked functions, we mark the (active) cells of
%  its same level within the support, and such that they are not contained
%  in the support of any function that has not been marked.
%
% Copyright (C) 2015, 2017 Eduardo M. Garau, Rafael Vazquez
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

function [ME, ind] = compute_cells_to_coarsen (hspace, hmsh, MF)

ME = cell (hmsh.nlevels, 1);
ind = cell (hmsh.nlevels, 1);

for lev = 1:hspace.nlevels
  if (~isempty(MF{lev}))
    ME{lev} = sp_get_cells (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), MF{lev});
    [ME{lev}, ~, ind{lev}] = intersect (ME{lev}, hmsh.active{lev});
    [funs, funs_per_cell] = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), ME{lev});
    funs_to_remain = setdiff (hspace.active{lev}, MF{lev});
    cells_to_remain = cellfun (@(x) any (ismember (x, funs_to_remain)), funs_per_cell);

    ME{lev} = ME{lev}(~cells_to_remain);
  end
end

end