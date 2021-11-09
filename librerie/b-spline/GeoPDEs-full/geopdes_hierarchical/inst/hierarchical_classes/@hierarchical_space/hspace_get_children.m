% HSPACE_GET_CHILDREN: compute the children of a given set of functions of the same level.
%
%     [children, flag, children_of_function] = hspace_get_children (hspace, lev, ind)
%
% Get the children of a given set of basis functions of level lev, with the
%  subdivision given by the "Proj" matrices
% All the children functions are stored in the same array.
%
% INPUT:
%
%     hspace: the hierarchical space (see hierarchical_space)
%     lev:    level of the functions to refine
%     ind:    indices of the functions in the tensor product space of level lev
%
% OUTPUT:
%
%     children: indices of the children, with the numbering of the tensor product space
%     flag:     a flag to tell whether all the input functions are active (1) 
%               active or deactivated (2), or if there is any passive function (0)
%     children_of_function: cell array with the children of each function
%
% Copyright (C) 2015, 2016, 2017 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2017-2019 Rafael Vazquez
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

function [children, flag, children_of_function] = hspace_get_children (hspace, lev, ind)

% % This computation would not work for the truncated basis
% Cmat = matrix_basis_change__ (hspace, lev+1);  
% [children,~] = find (Cmat(:,ind));
% children = unique (children);

if (isa (hspace.space_of_level(1), 'sp_scalar'))
  is_scalar = true;
elseif (isa (hspace.space_of_level(1), 'sp_vector'))
  is_scalar = false;
else
  error ('Unknown space type')
end
ndim = size (hspace.Proj{1}, 2);

z = cell (ndim, 1);
ind_sub = cell (ndim, 1);
children = [];
aux = cell (ndim, 1);
children_of_function = cell (numel(ind), 1);
if (is_scalar)
  [ind_sub{:}] = ind2sub ([hspace.space_of_level(lev).ndof_dir, 1], ind); % The extra 1 makes it work in any dimension

  for ii = 1:numel(ind_sub{1})
    for idim = 1:ndim
      aux{idim} = find (hspace.Proj{lev}{idim}(:,ind_sub{idim}(ii)));
    end
    [z{1:ndim}] = ndgrid (aux{:});
    auxI = sub2ind ([hspace.space_of_level(lev+1).ndof_dir, 1], z{:});
    children = union (children, auxI(:));
    children_of_function{ii} = auxI(:)';
  end
  
else
  cumsum_ndof_coarse = hspace.space_of_level(lev).cumsum_ndof;
  cumsum_ndof_fine = hspace.space_of_level(lev+1).cumsum_ndof;
  for icomp = 1:hspace.space_of_level(lev).ncomp_param
    ind_comp = ind(ind>cumsum_ndof_coarse(icomp) & ind<=cumsum_ndof_coarse(icomp+1)) - cumsum_ndof_coarse(icomp);
    [ind_sub{:}] = ind2sub ([hspace.space_of_level(lev).ndof_dir(icomp,:), 1], ind_comp); % The extra 1 makes it work in any dimension

    for ii = 1:numel(ind_sub{1})
      for idim = 1:ndim
        aux{idim} = find (hspace.Proj{lev}{icomp, idim}(:,ind_sub{idim}(ii)));
      end
      [z{1:ndim}] = ndgrid (aux{:});
      auxI = sub2ind ([hspace.space_of_level(lev+1).ndof_dir(icomp,:), 1], z{:});
      children = union (children, auxI(:)+cumsum_ndof_fine(icomp));
      children_of_function{ii} = auxI(:)'+cumsum_ndof_fine(icomp);
    end
  end  
end

if (nargout >= 2)
  flag = all (ismember (ind, hspace.active{lev}));
  if (~flag)
    flag = 2 *all (ismember (ind, union (hspace.active{lev}, hspace.deactivated{lev})));
  end
end

end
