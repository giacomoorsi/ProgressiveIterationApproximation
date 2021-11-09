% HMSH_REFINED_MESH_FOR_INTERFACE: generate an auxiliary (refined) hierarchical mesh, 
%  such that two adjacent elements on the interface are active on both patches
%
% USAGE:
%
%    [hmsh_aux, adjacent_elements, adjacent_elem_to_int_elem] = hmsh_refined_mesh_for_interface (hmsh, interface)
%
% INPUT:
%
%    hmsh:      a hierarchical multipatch mesh object (see hierarchical_mesh_mp)
%    interface: the interface information, for a single interface (as in hspace.interfaces(iref))
%
% OUTPUT:
%
%    hmsh_aux:           refined hierarchical mesh
%    adjacent_elements: for each level, and for each side of the interface,
%           indices of the adjacent active elements in hmsh_aux
%    adjacent_elem_to_int_elem: array of size (2, nelem_int). For each element on the interface (in hmsh_aux),
%           and for each side, global indices of the adjacent active elements (in hmsh)
%
% Copyright (C) 2017, 2018 Rafael Vazquez
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

function [hmsh_aux, adjacent_elements, adjacent_elem_to_int_elem] = hmsh_refined_mesh_for_interface (hmsh, interface)

  patch(1) = interface.patch1;
  patch(2) = interface.patch2;
  side(1) = interface.side1;
  side(2) = interface.side2;
  
  iedge = 0;
  hmsh_aux = hmsh;
  adjacent_elements = cell (hmsh.nlevels, 1);
  
  for lev = 1:hmsh.nlevels
    marked = cell (hmsh.nlevels, 1);
    Nelem = cumsum ([0, hmsh_aux.mesh_of_level(lev).nel_per_patch]);
    for ii = 1:2
      msh_patch_lev = hmsh_aux.mesh_of_level(lev).msh_patch{patch(ii)};
      nel_dir = msh_patch_lev.nel_dir;
%    ind = [1 1 2 2 3 3] in 3D, ind = [1 1 2 2] in 2D
%    ind2 = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
      ind = ceil (side(ii)/2);
      ind2 = setdiff (1:hmsh.ndim, ind);
      subindices = arrayfun (@(x) 1:x, nel_dir, 'UniformOutput', false);
      if (mod (side(ii), 2) == 1)
        subindices{ind} = 1;
      else
        subindices{ind} = nel_dir(ind);
      end
      [subindices{:}] = ndgrid (subindices{:});
      elems{ii} = reshape (sub2ind ([nel_dir, 1], subindices{:}), [nel_dir(ind2), 1]);
    end
    elems{1} = elems{1}(:)';
    elems{2} = reorder_elements (elems{2}, interface, nel_dir(ind2));
    [active_elements{1}, pos1] = ismember (elems{1}+Nelem(patch(1)), hmsh_aux.active{lev});
    [active_elements{2}, pos2] = ismember (elems{2}+Nelem(patch(2)), hmsh_aux.active{lev});
    
    interface_indices = active_elements{1} & active_elements{2};
    adjacent_elements{lev}{1} = hmsh_aux.active{lev}(pos1(interface_indices));
    adjacent_elements{lev}{2} = hmsh_aux.active{lev}(pos2(interface_indices));

    indices1 = (active_elements{1} & ~active_elements{2});
    indices2 = (active_elements{2} & ~active_elements{1});
    indices = union (pos1(indices1), pos2(indices2));
    marked{lev} = hmsh_aux.active{lev}(indices);
    hmsh_aux = hmsh_refine (hmsh_aux, marked);

    iedge_aux = iedge;
    Nelem_level = cumsum ([0 hmsh.nel_per_level]);

    for ii = 1:2
      iedge = iedge_aux;
      int_elems = adjacent_elements{lev}{ii};
      for iel = 1:numel(int_elems)
        elem_lev = int_elems(iel);
        iedge = iedge + 1;
        flag = 0;
        levj = lev;
        while (~flag)
          [is_active, pos] = ismember (elem_lev, hmsh.active{levj});
          if (is_active)
            adjacent_elem_to_int_elem(ii,iedge) = Nelem_level(levj) + pos;
            flag = 1;
          else
            elem_lev = hmsh_get_parent (hmsh, levj, elem_lev);
            flag = 0;
            levj = levj-1;
          end
        end
      end
    end

  end  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function elem = reorder_elements (elem, interface, nel_dir)
% Reorder elements of adjacent patches, to get a corresponding numbering
  ndim = numel (nel_dir) + 1;
  if (ndim == 2)
    if (interface.ornt == -1)
      elem = fliplr (elem(:)');
    else
      elem = elem(:)';
    end
  elseif (ndim == 3)
    if (interface.flag == -1)
      elem = elem';
    end
    if (interface.ornt1 == -1)
      elem = flipud (elem);
    end
    if (interface.ornt2 == -1)
      elem = fliplr (elem);
    end
    elem = elem(:)';
  end
end

