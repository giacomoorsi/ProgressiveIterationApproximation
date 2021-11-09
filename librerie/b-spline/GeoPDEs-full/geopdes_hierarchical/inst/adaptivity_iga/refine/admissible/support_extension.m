function [list_of_cells] = support_extension (hmsh, hspace, Q_ind, lev_Q, lev_s)

% SUPPORT_EXTENSION: compute the support extension of a given element, also with
%   respect to a different level
%
%   [list_of_cells] = support_extension (hmsh, hspace, Q_ind, lev_Q, lev_s)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the space of hierarchical splines (see hierarchical_space)
%   Q_ind:  indices of the input elements, all of the same level
%   lev_Q:  level of the elements in Q_ind
%   lev_s:  level for which we want to compute the support extension
%
% OUTPUT:
%
%   list_of_cells: indices, in the tensor product space, of the elements in 
%                   the support extension
%
% Copyright (C) 2017, 2018 Cesare Bracco, Rafael Vazquez
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

if (lev_Q == lev_s || nargin == 4)
    funs = sp_get_basis_functions (hspace.space_of_level(lev_Q), hmsh.mesh_of_level(lev_Q), Q_ind);
    list_of_cells = sp_get_cells (hspace.space_of_level(lev_Q), hmsh.mesh_of_level(lev_Q), funs);
elseif (lev_Q > lev_s)
    ancestors = hmsh_get_ancestors (hmsh, Q_ind, lev_Q, lev_s);
    list_of_cells = support_extension (hmsh, hspace, ancestors, lev_s, lev_s);
elseif (lev_Q < lev_s)
    descendants = hmsh_get_descendants (hmsh, Q_ind, lev_Q, lev_s);
    list_of_cells = support_extension (hmsh, hspace, descendants, lev_s, lev_s);
end

end