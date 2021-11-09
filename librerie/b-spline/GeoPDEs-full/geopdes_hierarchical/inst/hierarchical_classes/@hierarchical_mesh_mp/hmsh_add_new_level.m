% HMSH_ADD_NEW_LEVEL: add an empty level to the hierarchical mesh.
%
%   hmsh = hmsh_add_new_level (hmsh)
%
% INPUT:
%
%   hmsh:    object representing the hierarchical mesh (see hierarchical_mesh_mp)
%
% OUTPUT:
%
%   hmsh:   the object of the hierarchical mesh with one more level, which is empty
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

function hmsh = hmsh_add_new_level (hmsh)

  hmsh.nlevels = hmsh.nlevels + 1;
  hmsh.active{hmsh.nlevels} = [];
  hmsh.deactivated{hmsh.nlevels} = [];
  hmsh.nel_per_level(hmsh.nlevels) = 0;
  hmsh.mesh_of_level(hmsh.nlevels) = msh_refine (hmsh.mesh_of_level(hmsh.nlevels-1), hmsh.nsub);
  hmsh.msh_lev{hmsh.nlevels} = [];

end
