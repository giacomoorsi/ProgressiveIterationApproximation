% ADAPTIVITY_MARK_COARSENING: mark cells or basis functions for coarsening according to 
%  the marking strategy in adaptivity_data.mark_strategy, taking into account 
%  the computed error estimators in the variable est.
%
%   [marked, nmarked] = adaptivity_mark_coarsening (est, hmsh, hspace, adaptivity_data)
%
% INPUT:
%
%   est:     result obtained from the estimate step (see adaptivity_estimate_laplace)
%   hmsh:    object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace:  object representing the coarse space of hierarchical splines (see hierarchical_space)
%   adaptivity_data: a structure with the data for the adaptivity method.
%                    In particular, it must contain the following fields:
%      -mark_strategy: the possible marking strategies are:
%      -mark_param:    parameter for marking, 0 < mark_param < 1.
%      -flag:          elements or functions, according to est
%
% OUTPUT:
%
%    marked:  cell-array with the indices of marked cells (or functions) for each level, in the tensor-product setting
%    nmarked: total number of marked cells (or functions)
%
% Copyright (C) 2015, 2016, 2017 Eduardo M. Garau, Rafael Vazquez
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


function [marked, nmarked] = adaptivity_mark_coarsening (est, hmsh, hspace, adaptivity_data)

switch adaptivity_data.flag
    case 'elements' %, disp ('marking elements for refinement')
    case 'functions' %, disp ('marking basis functions for refinement')
    otherwise, error ('adaptivity_mark: Unknown option %s', adaptivity_data.flag)
end

max_est = max (est);
aux_marked = zeros (size (est));

switch adaptivity_data.mark_strategy
 case 'GR'
  aux_marked = ones (size (est));
 case {'MS', 'GERS'}
  nn = round(adaptivity_data.mark_param_coarsening*numel(est));
  [~, ind] = sort(est);
  aux_marked(ind(1:nn)) = 1;
% case 'MS'
%  aux_marked(est < adaptivity_data.mark_param_coarsening * max_est) = 1;
% case 'GERS'
%  est_sum2 = sum (est.^2);
%  [est2_ordered, perm] = sort (est.^2, 'ascend');
%  index = find (cumsum (est2_ordered.^2) < adaptivity_data.mark_param^2 * est_sum2, 1, 'last');
% %  index = find (est2_ordered < 1.001 * est2_ordered(index), 1, 'first');
%  aux_marked(perm(1:index)) = 1;
end

marked_list = find (aux_marked);
nmarked = numel (marked_list);

marked = cell (hmsh.nlevels, 1);

switch (lower (adaptivity_data.flag))
  case 'elements'
    aux = cumsum ([0, hmsh.nel_per_level]);
    for lev = 1:hmsh.nlevels
      elems = aux(lev)+1:aux(lev+1);
      [~,ind,~] = intersect (elems, marked_list);
      marked{lev} = hmsh.active{lev}(ind);
    end
  case 'functions'
    aux = cumsum ([0, hspace.ndof_per_level]);
    for lev = 1:hmsh.nlevels
      funs = aux(lev)+1:aux(lev+1);
      [~,ind,~] = intersect (funs, marked_list);
      marked{lev} = hspace.active{lev}(ind);
    end
end

end
