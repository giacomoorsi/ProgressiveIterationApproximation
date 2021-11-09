% HSPACE_REMOVE_EMPTY_LEVELS: remove the last levels of the hierarchical space, if empty.
%   The finest level will be the same as in the hierarchical mesh object,
%   even if there are no active functions.
%
%   hspace = hmsh_remove_empty_levels (hspace, hmsh)
%
% INPUT:
%
%   hspace: object representing the hierarchical space (see hierarchical_space_mp)
%   hmsh:  object representing the hierarchical mesh, with the last level already removed (see hierarchical_mesh_mp)
%
% OUTPUT:
%
%   hspace: the object of the hierarchical space after removing the empty levels
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

function hspace = hspace_remove_empty_levels (hspace, hmsh)

if (~isempty (hmsh.boundary))
  for iside = 1:numel(hmsh.boundary)
    hmsh.boundary(iside) = hmsh_remove_empty_levels (hmsh.boundary(iside));
    hspace.boundary(iside) = hspace_remove_empty_levels (hspace.boundary(iside), hmsh.boundary(iside));
  end
end

for ilev = hspace.nlevels:-1:hmsh.nlevels+1
  if (isempty (hspace.active{ilev}))
    hspace.space_of_level(ilev) = [];
    hspace.Proj(ilev-1,:) = [];

    hspace.active(ilev) = [];
    hspace.deactivated(ilev) = [];
    hspace.ndof_per_level(ilev) = [];
    hspace.nlevels = hspace.nlevels - 1;
  end
end

end
