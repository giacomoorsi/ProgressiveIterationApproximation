% HMSH_REFINE: refine the hierarchical mesh, updating the fields of the object.
%
%   [hmsh, new_elements] = hmsh_refine (hmsh, marked)
%
% INPUT:
%
%   hmsh:    object representing the coarse hierarchical mesh (see hierarchical_mesh_mp)
%   marked:  cell array with the indices, in the Cartesian grid, of the marked elements for each level
%
% OUTPUT:
%
%   hmsh:         object representing the refined hierarchical mesh (see hierarchical_mesh_mp)
%   new_elements: cell array with the global indices of the new active elements for each level
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function [hmsh, new_elements] = hmsh_refine (hmsh, M)

boundary = ~isempty (hmsh.boundary);

% Computation of a new grid if a new level is activated
if (~isempty(M{hmsh.nlevels}))
  hmsh = hmsh_add_new_level (hmsh);
end

% Update the set of active elements
old_elements = hmsh.active;
[hmsh, new_elements] = update_active_cells (hmsh, M);

% Update msh_lev
hmsh.msh_lev = update_msh_lev (hmsh, old_elements, new_elements);
% hmsh.msh_lev = cell (hmsh.nlevels,1);
% for ilev = 1 : hmsh.nlevels
%   hmsh.msh_lev{ilev} = msh_evaluate_element_list (hmsh.mesh_of_level(ilev), hmsh.active{ilev});
% end

% Update the boundary, calling the function recursively
if (boundary)
  if (hmsh.ndim > 1)
    M_boundary = cell (size (M));
    for lev = 1:min(numel (M), hmsh.boundary.nlevels)
      Nelem = cumsum ([0 hmsh.mesh_of_level(lev).nel_per_patch]);
      Nelem_bnd = cumsum ([0 hmsh.boundary.mesh_of_level(lev).nel_per_patch]);
      for iptc = 1:hmsh.boundary.npatch
        patch_number = hmsh.boundary.mesh_of_level(1).patch_numbers(iptc);
        side_number  = hmsh.boundary.mesh_of_level(1).side_numbers(iptc);
        [~,indices,~] = intersect (Nelem(patch_number)+1:Nelem(patch_number+1), M{lev});
        bnd_indices = get_boundary_indices (side_number, hmsh.mesh_of_level(lev).msh_patch{patch_number}.nel_dir, indices);
        M_boundary{lev} = union (M_boundary{lev}, bnd_indices+Nelem_bnd(iptc));
      end
    end
    hmsh.boundary = hmsh_refine (hmsh.boundary, M_boundary);
  end
else
  hmsh.boundary = [];
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hmsh, new_cells] = update_active_cells (hmsh, M)
%
% function [hmsh, new_cells] = update_active_cells (hmsh, M)
%
% Update the sets of active cells (hmsh.active) and deactivated cells (hmsh.deactivated) in each level when
% refining the cells in M. This function also updates hmsh.nlevels, hmsh.nel and hmsh.nel_per_level
%
% INPUT
%     hmsh: object representing the coarse hierarchical mesh (see hierarchical_mesh_mp)
%     M{lev}: global indices of marked cells of level lev (one row per cell)
%
% OUTPUT
%     hmsh: object representing the refined hierarchical mesh (see hierarchical_mesh_mp)
%     new_cells{lev}: global indices of the new cells of level lev (one row per cell)
%

% nlevels = hmsh.nlevels;
nlevels = numel (M);

new_cells = cell (hmsh.nlevels, 1);

% Deactivate the cells to be refined, and compute their children
for lev = 1:nlevels
  if (~isempty(M{lev}))
    [dummy, indE] = ismember (M{lev}, hmsh.active{lev});
%       if (~all (dummy))
%         warning('update_active_cells: Some nonactive cells were selected');
%       end
    hmsh.active{lev}(indE) = [];
    hmsh.deactivated{lev} = union (hmsh.deactivated{lev}, M{lev});
      
    new_cells{lev+1} = hmsh_get_children (hmsh, lev, M{lev});
    hmsh.active{lev+1} = union (hmsh.active{lev+1}, new_cells{lev+1});
  end
end

% Update hmsh.nel_per_level and hmsh.nel
hmsh.nel_per_level = cellfun (@numel, hmsh.active);
hmsh.nel = sum (hmsh.nel_per_level);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msh_lev = update_msh_lev (hmsh, old_elements, new_elements)
%
% function msh_lev = update_msh_lev (hmsh, old_elements, new_elements)
%
% Update the information in msh_lev, computing only the elements that have been added to the mesh
%
% INPUT
%     hmsh: object representing the fine hierarchical mesh (see hierarchical_mesh_mp)
%     old_elements{lev}: active elements in the previous coarse mesh
%     new_elements{lev}: elements that have been added after refinement
%
% OUTPUT
%     msh_lev: the structures with the msh information for each level
%

msh_lev = cell (hmsh.nlevels, 1);

for lev = 1:hmsh.nlevels
  if (lev > numel (old_elements) || numel (old_elements{lev}) == 0)
    msh_lev{lev} = msh_evaluate_element_list (hmsh.mesh_of_level(lev), hmsh.active{lev});
  else
    [~, iold, iold_act] = intersect (old_elements{lev}, hmsh.active{lev});
    msh_lev{lev}.npatch = hmsh.npatch;
    msh_lev{lev}.ndim = hmsh.ndim;
    msh_lev{lev}.rdim = hmsh.rdim;
    msh_lev{lev}.nel = hmsh.nel_per_level(lev);
    msh_lev{lev}.elem_list = hmsh.active{lev}(:).';

    if (isempty (new_elements{lev}))
      indices = iold_act;
      msh_new = struct ('quad_weights', [], 'geo_map', [], 'geo_map_jac', [], 'geo_map_der2', [], 'jacdet', [], 'element_size', []);
    else
      msh_new = msh_evaluate_element_list (hmsh.mesh_of_level(lev), new_elements{lev});
      [~, ~, inew_act] = intersect (new_elements{lev}, hmsh.active{lev});
      indices = [iold_act(:); inew_act(:)];
    end
    msh_lev{lev}.quad_weights(:,indices) = [hmsh.msh_lev{lev}.quad_weights(:,iold), msh_new.quad_weights];
    msh_lev{lev}.geo_map(:,:,indices) = cat (3, hmsh.msh_lev{lev}.geo_map(:,:,iold), msh_new.geo_map);
    msh_lev{lev}.geo_map_jac(:,:,:,indices) = cat (4, hmsh.msh_lev{lev}.geo_map_jac(:,:,:,iold), msh_new.geo_map_jac);
    msh_lev{lev}.geo_map_der2(:,:,:,:,indices) = cat (5, hmsh.msh_lev{lev}.geo_map_der2(:,:,:,:,iold), msh_new.geo_map_der2);
    msh_lev{lev}.jacdet(:,indices) = [hmsh.msh_lev{lev}.jacdet(:,iold), msh_new.jacdet];
    msh_lev{lev}.element_size(:,indices) = [hmsh.msh_lev{lev}.element_size(:,iold), msh_new.element_size];

    Nelem = cumsum ([0 hmsh.mesh_of_level(lev).nel_per_patch]);
    for iptc = 1:msh_lev{lev}.npatch
      [~,indices,~] = intersect (Nelem(iptc)+1:Nelem(iptc+1), msh_lev{lev}.elem_list);
      msh_lev{lev}.nel_per_patch(iptc) = numel (indices);
      msh_lev{lev}.elem_list_of_patch{iptc} = indices;
      msh_lev{lev}.nel_dir_of_patch{iptc} = hmsh.mesh_of_level(lev).msh_patch{iptc}.nel_dir;
    end
    msh_lev{lev}.nqn = hmsh.mesh_of_level(lev).msh_patch{1}.nqn;
    msh_lev{lev}.nqn_dir = hmsh.mesh_of_level(lev).msh_patch{1}.nqn_dir;
  end
end

end
