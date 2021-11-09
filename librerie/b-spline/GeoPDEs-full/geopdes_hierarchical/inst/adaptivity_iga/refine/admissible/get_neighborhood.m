function [ list_of_cells ] = get_neighborhood (hmsh, hspace, Q_ind, lev_Q, m, adm_type)

% GET_NEIGHBORHOOD: compute the neighborhood of the cells of a hierarchical mesh
%
%   [list_of_cells] = get_neighborhood (hmsh, hspace, Q_ind, Q_lev, m, adm_type)
%
% INPUT:
%
%   hmsh:     object representing the hierarchical mesh (see hierarchical_mesh)
%   hspace:   object representing the space of hierarchical splines (see hierarchical_space)
%   Q_ind:    indices of the input elements, all of the same level
%   Q_lev:    level of the elements in Q_ind
%   m:        admissibility class of the hierarchical mesh, an integer number
%   adm_type: admissibility type, either 'T-admissible' or 'H-admissible'
%
% OUTPUT:
%
%   list_of_cells: cell array with the indices, in the tensor product space,
%      of the active elements in the neighborhood of the elements in Q_ind
%           
%  The definition of the neighborhood is from the papers: 
%      A. Buffa and C. Giannelli
%      Adaptive isogeometric methods with hierarchical splines: error
%       estimator and convergence
%      Math. Models Meth. Appl. Sci., 2016
%
%      C. Bracco, C. Giannelli, R. Vazquez
%      Refinement algorithms for adaptive isogeometric methods with
%      hierarchical splines, Axioms, 2018
%
%
% Copyright (C) 2017, 2018, 2019 Cesare Bracco, Rafael Vazquez
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

if (strcmpi (adm_type, 'T-admissible'))
  lev_s = lev_Q-m+2; 
elseif (strcmpi (adm_type, 'H-admissible'))
  lev_s = lev_Q-m+1;
else
  error ('get_neighborhood: unknown type of admissibility')
end
lev_n=lev_Q-m+1;  %level of the elements in the neighborhood

if (lev_n < 1)
  list_of_cells = [];
else
  supp_ext = support_extension (hmsh, hspace, Q_ind, lev_Q, lev_s);
  if (strcmpi (adm_type, 'T-admissible'))
    ancestors = hmsh_get_parent (hmsh, lev_s, supp_ext);
    list_of_cells = intersect (ancestors, hmsh.active{lev_n});
  else
    list_of_cells = intersect (supp_ext, hmsh.active{lev_n});
  end
end

end