% HMSH_GET_PARENT: compute the parent cell of a given cell of a certain level.
%
%     [parents, flag] = hmsh_get_parent (hmsh, lev, ind)
%
% Get the parent cell of a given cell of level lev, with the subdivision given by hmsh.nsub.
%  If several cells are given as input, return as output all their parent
%  cells in a unique array, without repetitions.
%
% INPUT:
%
%     hmsh: the hierarchical mesh (see hierarchical_mesh_mp)
%     lev:  level of the subdivided cells to compute their parents
%     ind:  indices of the cells in the global multipatch grid of that level
%
% OUTPUT:
%
%     parent: index of the parent, with the numbering of the multipatch grid
%     flag:   a flag to tell whether all the input cells are active (1) 
%               active or deactivated (2), or if there is any passive cell (0)
%
% Copyright (C) 2016 Eduardo M. Garau, Rafael Vazquez
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

function [parent, flag] = hmsh_get_parent (hmsh, lev, ind)

if (lev < 2 || lev > hmsh.nlevels)
  error ('The level should be between 2 and the number of levels of the mesh')
end
if (any (ind > hmsh.mesh_of_level(lev).nel))
  error ('There are some indices greater than the number of elements of the level')
end

parent = [];
ndim = hmsh.ndim;
nsub = hmsh.nsub;

Nelem = cumsum ([0 hmsh.mesh_of_level(lev).nel_per_patch]);
Nelem_coarse = cumsum ([0 hmsh.mesh_of_level(lev-1).nel_per_patch]);
aux = cell (ndim, 1);
for iptc = 1:hmsh.npatch
  [~,indices,~] = intersect (Nelem(iptc)+1:Nelem(iptc+1), ind);
  z = cell (ndim, 1);
  cells_sub = cell (ndim, 1);
  [cells_sub{:}] = ind2sub ([hmsh.mesh_of_level(lev).msh_patch{iptc}.nel_dir, 1], indices); % The extra 1 makes it work in any dimension

  for ii = 1:numel(cells_sub{1})
    for idim = 1:ndim
      aux{idim} = floor ((cells_sub{idim}(ii) + nsub(idim) - 1) / nsub(idim));
    end
    [z{1:ndim}] = ndgrid (aux{:});
    auxI = sub2ind ([hmsh.mesh_of_level(lev-1).msh_patch{iptc}.nel_dir, 1], z{:});
    parent = union (parent, auxI(:)+Nelem_coarse(iptc));
  end
end

if (nargout == 2)
  flag = all (ismember (ind, hmsh.active{lev}));
  if (~flag)
    flag = 2 *all (ismember (ind, union (hmsh.active{lev}, hmsh.deactivated{lev})));
  end
end

end
