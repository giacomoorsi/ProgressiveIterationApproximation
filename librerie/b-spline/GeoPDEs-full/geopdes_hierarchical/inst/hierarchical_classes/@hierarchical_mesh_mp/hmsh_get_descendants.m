function [descendants, descendants_of_cell] = hmsh_get_descendants (hmsh, Q_ind, lev_Q, lev)

% HMSH_GET_DESCENDANTS: compute the descendants of a given list of elements of a hierarchical mesh
%
%   [descendants, descendants_of_cell] = hmsh_get_descendants (hmsh, Q_ind, lev_Q, lev)
%
% INPUT:
%
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%   Q_ind:  indices of the input elements, all of the same level
%   Q_lev:  level of the elements in Q_ind
%   lev:    level for which to compute the ancestors (lev_Q < lev)
%
% OUTPUT:
%
%   descendants: indices of the descendants of level lev
%   descendants_of_cell: indices of the descendants (rows) for each cell in ind (column)
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

if (lev_Q < 1 || lev_Q > hmsh.nlevels-1)
  error ('The level should be between 1 and the number of levels of the mesh minus one')
end

if (lev <= lev_Q)
  error('lev<= lev_Q: use get_parent')
end

if (lev > hmsh.nlevels)
  error ('The level of the descendants is higher than the finest level of the mesh')
end
  
if (nargout == 2)
  if (lev == lev_Q+1)
    [descendants,~,descendants_of_cell] = hmsh_get_children (hmsh, lev_Q, Q_ind);
  else
    [~,~,children_of_cell] = hmsh_get_children (hmsh, lev_Q, Q_ind);
    [descendants, descendants_of_children] = hmsh_get_descendants (hmsh, children_of_cell(:), lev_Q+1, lev);
    descendants_of_cell = reshape (descendants_of_children, [], numel(Q_ind));
  end  
else
  if (lev == lev_Q+1)
    descendants = hmsh_get_children (hmsh, lev_Q, Q_ind);
  else
    children = hmsh_get_children (hmsh, lev_Q, Q_ind);
    descendants = hmsh_get_descendants (hmsh, children, lev_Q+1, lev);
  end
end

end
