function genealogy = hmsh_get_ancestors (hmsh, Q_ind, lev_Q, lev)

% HMSH_GET_ANCESTORS: compute the ancestors of a given list of elements of a hierarchical mesh
%
%   ancestors = hmsh_get_ancestors (hmsh, Q_ind, lev_Q, lev)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   Q_ind:  indices of the input elements, all of the same level
%   Q_lev:  level of the elements in Q_ind
%   lev:    level for which to compute the ancestors (lev_Q > lev)
%
% OUTPUT:
%
%   ancestors: indices of the ancestors of level lev
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


if lev >= lev_Q
    error('lev>= lev_Q: use get_children')
end
  
if lev == lev_Q-1
    genealogy = hmsh_get_parent (hmsh, lev_Q, Q_ind);
else
    genealogy = hmsh_get_ancestors (hmsh, hmsh_get_parent(hmsh, lev_Q, Q_ind), lev_Q-1, lev);
end
end