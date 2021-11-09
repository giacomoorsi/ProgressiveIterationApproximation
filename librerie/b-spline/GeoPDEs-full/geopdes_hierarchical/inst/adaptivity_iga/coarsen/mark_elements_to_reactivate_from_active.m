% MARK_ELEMENTS_TO_REACTIVATE_FROM_ACTIVE: compute the elements to
%     reactivate from the set of marked elements
%
%   [deact_marked, num] = mark_elements_to_reactivate_from_active (marked, hmsh, hspace, adaptivity_data)
%
% INPUT:
%
%   marked:  cell-array with the indices of marked cells (or functions) for each level, in the tensor-product setting
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%   adaptivity_data: a structure with the data for the adaptivity method. In particular, it contains the fields:
%     - coarsening_flag: either 'all' (default) or 'any', to decide how
%                        many children must be marked to reactivate an element.
%     - adm_class: admissibility class, an integer value, zero by default (no admissibility).
%     - adm_type:  admissibility_type, either 'T-admissible' (default) or 'H-admissible'.
%
% OUTPUT:
%
%    deact_marked:  cell-array with the indices of marked cells (or functions) for each level to be reactivated, in the tensor-product setting
%    num         :  number of cells (or elements) to be reactivated
%
%    The coarsening algorithm is detailed in the paper
%      M. Carraturo, C. Giannelli, A. Reali, R. Vazquez
%      Suitably graded THB-spline refinement and coarsening: Towards 
%      an adaptive isogeometric analysis of additive manufacturing processes. 
%      Comput. Methods Appl. Mech. Engrg., 2019.
%
% Copyright (C) 2016, 2017, 2018 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2018, 2019 Massimo Carraturo, Rafael Vazquez
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

function [deact_marked, num] = mark_elements_to_reactivate_from_active (marked, hmsh, hspace, adaptivity_data)

if (~isfield (adaptivity_data, 'coarsening_flag'))
  adaptivity_data.coarsening_flag = 'all';
end

if (~isfield(adaptivity_data,'adm_class') || adaptivity_data.adm_class < 2)
  adm_class = 0;
else
  adm_class = adaptivity_data.adm_class;
end

if (~isfield(adaptivity_data,'adm_type') || isempty (adaptivity_data.adm_type))
  adm_type = 'T-admissible';
else
  adm_type = adaptivity_data.adm_type;
end


deact_marked = cell (hmsh.nlevels, 1);

for lev = hmsh.nlevels-1:-1:1
  if (~isempty(marked{lev+1}))
    [parents, flag] = hmsh_get_parent (hmsh, lev+1, marked{lev+1});
    if (flag ~= 1)
      error ('Some nonactive elements were marked.')
    end

    [~,~,children_per_cell] = hmsh_get_children (hmsh, lev, parents);    

    ind = all (ismember (children_per_cell, hmsh.active{lev+1}));
    parents = parents(ind);
    children_per_cell = children_per_cell(:,ind);

    if (strcmpi (adaptivity_data.coarsening_flag, 'any'))
      deact_marked{lev} = parents;
    elseif (strcmpi (adaptivity_data.coarsening_flag, 'all'))
      ind2 = all (ismember (children_per_cell, marked{lev+1}));
      parents = parents(ind2);
      children_per_cell = children_per_cell(:,ind2);
      deact_marked{lev} = parents;
    else
      error ('Unknown option for coarsening, in adaptivity_data.coarsening_flag')
    end
    marked{lev+1} = children_per_cell;
    
% Algorithm to recover admissible meshes
    if (adm_class)
      lev_s = lev + adm_class;
      if (lev_s > hmsh.nlevels)
        continue
      else
        active_and_deact = union (hmsh.active{lev_s}, hmsh.deactivated{lev_s});
        active_and_deact = setdiff (active_and_deact, marked{lev_s});
        keep_inds = [];
        
        if (strcmpi (adm_type, 'T-admissible'))
          supp_ext = support_extension (hmsh, hspace, children_per_cell(:), lev+1, lev+1);
          [~, descendants_of_cell] = hmsh_get_descendants (hmsh, supp_ext, lev+1, lev_s);
          for iel = 1:numel(deact_marked{lev})
            supp_ext_local = support_extension (hmsh, hspace, children_per_cell(:,iel), lev+1, lev+1);
            [~,ia,~] = intersect (supp_ext, supp_ext_local);
            if (isempty (intersect (descendants_of_cell(:,ia), active_and_deact)))
              keep_inds = [keep_inds, iel];
            end          
          end
          
        elseif (strcmpi (adm_type, 'H-admissible'))
          supp_ext = support_extension (hmsh, hspace, deact_marked{lev}, lev, lev);
          [~, descendants_of_cell] = hmsh_get_descendants (hmsh, supp_ext, lev, lev_s);
          for iel = 1:numel(deact_marked{lev})
            supp_ext_local = support_extension (hmsh, hspace, deact_marked{lev}(iel), lev, lev);
            [~,ia,~] = intersect (supp_ext, supp_ext_local);
            if (isempty (intersect (descendants_of_cell(:,ia), active_and_deact)))
              keep_inds = [keep_inds, iel];
            end          
          end
        end
        deact_marked{lev} = deact_marked{lev}(keep_inds);
      end
    end
    
  end
end
  
num = sum (cellfun (@numel, deact_marked));
    
end
