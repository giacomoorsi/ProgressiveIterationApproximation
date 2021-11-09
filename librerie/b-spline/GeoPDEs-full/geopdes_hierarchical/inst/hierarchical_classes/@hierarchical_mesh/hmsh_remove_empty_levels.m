% HMSH_REMOVE_EMPTY_LEVELS: remove the last levels of the hierarchical mesh, only if empty.
%
%   hmsh = hmsh_remove_empty_levels (hmsh)
%
% INPUT:
%
%   hmsh:    object representing the hierarchical mesh (see hierarchical_mesh)
%
% OUTPUT:
%
%   hmsh:   the object of the hierarchical mesh after removing the empty levels.
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
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

function hmsh = hmsh_remove_empty_levels (hmsh)

if (~isempty (hmsh.boundary))
  if (hmsh.ndim > 1)
    for iside = 1:numel(hmsh.boundary)
      hmsh.boundary(iside) = hmsh_remove_empty_levels (hmsh.boundary(iside));
    end
  end
end

for ilev = hmsh.nlevels:-1:2
  if (isempty (hmsh.active{ilev}))
    hmsh.nlevels = hmsh.nlevels - 1;
    hmsh.active = hmsh.active(1:hmsh.nlevels);
    hmsh.deactivated = hmsh.deactivated(1:hmsh.nlevels);
    hmsh.nel_per_level = hmsh.nel_per_level(1:hmsh.nlevels);
    hmsh.mesh_of_level = hmsh.mesh_of_level(1:hmsh.nlevels);
    hmsh.msh_lev = hmsh.msh_lev(1:hmsh.nlevels);    
  else
    break
  end
end

end
