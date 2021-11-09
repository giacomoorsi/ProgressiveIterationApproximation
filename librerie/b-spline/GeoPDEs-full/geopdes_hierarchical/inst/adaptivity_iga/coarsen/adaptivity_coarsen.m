% ADAPTIVITY_COARSEN: coarsen the hierarchical mesh and space, updating the corresponding structures hmsh and hspace.
%  The refinement can be done marking either elements or basis functions.
%
%   [hmsh, hspace] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data)
%
% INPUT:
%
%   hmsh:   object representing the fine hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the fine space of hierarchical splines (see hierarchical_space)
%   marked: cell array with the indices, in the tensor product space, of the marked elements/functions
%            to reactivate for each level
%   adaptivity_data: a structure with the data for the adaptivity method.
%                    In particular, it contains the field 'flag', that can take the value
%                    'elements' or 'functions', depending on the coarsening strategy.
%
% OUTPUT:
%
%   hmsh:   object representing the coarsened hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarsened space of hierarchical splines (see hierarchical_space)
%
% Copyright (C) 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2018-2019 Rafael Vazquez
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

function [hmsh, hspace, Ccoar] = adaptivity_coarsen (hmsh, hspace, marked, adaptivity_data)

switch (adaptivity_data.flag)
  case 'functions'
    error ('Coarsening by functions is not implemented yet')
  case 'elements'
    marked_elements = marked;
end
[reactivated_elements, ~] = mark_elements_to_reactivate_from_active (marked_elements, hmsh, hspace, adaptivity_data);

hmsh_fine = hmsh;
[hmsh, removed_cells] = hmsh_coarsen (hmsh, reactivated_elements);

reactivated_fun = functions_to_reactivate_from_cells (hmsh, hspace, reactivated_elements);

if (nargout == 3)
  hspace_fine = hspace;
  hspace = hspace_coarsen (hspace, hmsh, reactivated_fun, removed_cells);
  warning ('Coarsening matrix computed with an expensive L^2 projection')
  M = op_u_v_hier (hspace, hspace, hmsh);
  G = op_u_v_hier (hspace_fine, hspace_in_finer_mesh(hspace, hmsh, hmsh_fine), hmsh_fine);
  Ccoar = M \ G; Ccoar(abs(Ccoar) < 1e-12) = 0;
else
  hspace = hspace_coarsen (hspace, hmsh, reactivated_fun, removed_cells);
end
  
hmsh = hmsh_remove_empty_levels (hmsh);
hspace = hspace_remove_empty_levels (hspace, hmsh);

end
