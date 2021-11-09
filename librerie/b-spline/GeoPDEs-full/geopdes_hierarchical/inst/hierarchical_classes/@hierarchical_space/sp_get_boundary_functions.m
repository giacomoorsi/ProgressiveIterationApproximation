% SP_GET_BOUNDARY_FUNCTIONS: indices of the degrees of freedom of the given boundaries.
%
%   dofs = sp_get_boundary_functions (hspace, [hmsh], sides)
%
% INPUT:
%
%  hspace: object representing the hierarchical space of trial functions (see hierarchical_space)
%  hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
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

dofs = [];
for iside = drchlt_sides
  dofs = union (dofs, hspace.boundary(iside).dofs);
end
    
end
