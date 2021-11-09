% CHANGE_CONNECTIVITY_LOCALIZED_CSUB: compute the new connectivity related to the localized Csub.
%
% function sp_lev = change_connectivity_localized_Csub (sp_lev, hspace, level)
%
% Compute the new indices and number of dofs related to the localized Csub
%
% INPUT:  
%   sp_lev: object representing the current space of level (see space_scalar / space_vector)
%   hspace: object representing the hierarchical space (see hierarchical_space)
%   level:  level of the input space of fixed level
%
% OUTPUT:
%   sp_lev: object representing the current space of level (see space_scalar / space_vector)
%            with a modified connectivity
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2017-2019 Rafael Vazquez
% Copyright (C) 2018-2019 Luca Coradello
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

function sp_lev = change_connectivity_localized_Csub (sp_lev, hspace, ilev)

row_indices = hspace.Csub_row_indices{ilev};
[~,position] = ismember (sp_lev.connectivity, row_indices);
if (any (~position))
  error ('The given indices do not contain all the functions.')
end
sp_lev.ndof = numel (row_indices);
sp_lev.connectivity = position;

end
