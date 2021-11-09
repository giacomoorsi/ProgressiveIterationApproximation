% ADAPTIVITY_REFINE: refine the hierarchical mesh and space, updating the corresponding structures hmsh and hspace.
%  The refinement can be done marking either elements or basis functions.
%
%   [hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data)
%
% INPUT:
%
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%   marked: cell array with the indices, in the tensor product space, of the marked elements/functions
%            for each level
%   adaptivity_data: a structure with the data for the adaptivity method.
%                    It contains the field 'flag', that can take the value
%                    'elements' or 'functions', depending on the refinement strategy.
%
% OUTPUT:
%
%   hmsh:   object representing the refined hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the refined space of hierarchical splines (see hierarchical_space)
%   Cref:   refinement matrix, to pass from the original space to the refined one
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

function [hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data)

switch (adaptivity_data.flag)
  case 'functions'
    marked_elements = compute_cells_to_refine (hspace, hmsh, marked);
  case 'elements'
    marked_elements = marked;
end

if (isfield (adaptivity_data, 'adm_class'))
  marked_elements = mark_admissible (hmsh, hspace, marked_elements, adaptivity_data);
  [hmsh, new_cells] = hmsh_refine (hmsh, marked_elements);
  marked_functions = compute_functions_to_deactivate (hmsh, hspace, marked_elements, 'elements');
else
  [hmsh, new_cells] = hmsh_refine (hmsh, marked_elements);
  marked_functions = compute_functions_to_deactivate (hmsh, hspace, marked, adaptivity_data.flag);  
end

if (nargout == 3)
  [hspace, Cref] = hspace_refine (hspace, hmsh, marked_functions, new_cells);
else
  hspace = hspace_refine (hspace, hmsh, marked_functions, new_cells);
end

end