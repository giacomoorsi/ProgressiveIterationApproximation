% HSPACE_ADMISSIBILITY_CLASS: compute the admissibility class of the hierarchical mesh.
%
%   [adm_class] = hspace_admissibility_class (hspace, hmsh, [adm_type])
%
% INPUT:
%
%   hspace:   object representing the coarse space of hierarchical splines (see hierarchical_space)
%   hmsh:     object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   adm_type: either 'T-admissible' or 'H-admissible'. If not given 
%              the function will use the admissibility type corresponding
%              to the chosen basis (truncated or not).
%
% OUTPUT:
%
%   adm_class: integer value with the admissibility class.
%
% Copyright (C) 2018-2019 Cesare Bracco, Rafael Vazquez
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
%returns 1 if the hierarchical mesh is admissible of class m, 0 otherwise


function adm_class = hspace_admissibility_class (hspace, hmsh, adm_type)

if (nargin == 3)
  if (strcmpi(adm_type, 'H-admissible'))
    recompute = hspace.truncated;    
  elseif (strcmpi(adm_type, 'T-admissible'))
    recompute = ~hspace.truncated;
  else
    error ('Unknown admissibility type')
  end
  if (recompute)
    hspace.truncated = ~hspace.truncated;
    [hspace.Csub, hspace.Csub_row_indices] = hspace_subdivision_matrix (hspace, hmsh);
  end
end

inds = cell (hmsh.nel, 1);
minimum_level = zeros (1, hmsh.nel);
maximum_level = zeros (1, hmsh.nel);
% element_level = zeros (1, hmsh.nel);
shifting_index = cumsum ([0 hmsh.nel_per_level]);

for lev = 1:hmsh.nlevels
%   element_level(shifting_index(lev)+1:shifting_index(lev+1)) = lev;
  [~,funs] = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), hmsh.active{lev});
  for iel = 1:hmsh.nel_per_level(lev)
    el_number = iel+shifting_index(lev);
    [~,position] = ismember (funs{iel}, hspace.Csub_row_indices{lev});
    [~,inds{el_number}] = find (hspace.Csub{lev}(position,:));
    minimum_level(el_number) = find (min (inds{el_number}) < cumsum(hspace.ndof_per_level), 1);
    maximum_level(el_number) = find (max (inds{el_number}) <= cumsum(hspace.ndof_per_level), 1);
  end
end

adm_class = max (maximum_level - minimum_level) + 1;

end
