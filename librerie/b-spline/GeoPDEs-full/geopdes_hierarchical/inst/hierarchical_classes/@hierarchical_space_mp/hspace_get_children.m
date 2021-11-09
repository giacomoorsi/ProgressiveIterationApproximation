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
%     hspace: the hierarchical space (see hierarchical_space_mp)
%     lev:    level of the functions to refine
%     ind:    indices of the functions in the multipatch space of level lev
%
% OUTPUT:
%
%     children: indices of the children, with the numbering of the tensor product space
%     flag:     a flag to tell whether all the input functions are active (1) 
%               active or deactivated (2), or if there is any passive function (0)
%     children_of_function: cell array with the children of each function
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

function [children, flag, children_of_function] = hspace_get_children (hspace, lev, ind)

% % This computation would not work for the truncated basis
% Cmat = matrix_basis_change__ (hspace, lev+1);  
% [children,~] = find (Cmat(:,ind));
% children = unique (children);

if (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_scalar'))
  is_scalar = true;
  ndim = size (hspace.Proj{1}, 2);
elseif (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_vector'))
  is_scalar = false;
  ndim = size (hspace.Proj{1}, 2);
else
  error ('Unknown space type')
end

z = cell (ndim, 1);
ind_sub = cell (ndim, 1);
children = [];

npatch = numel (hspace.space_of_level(1).sp_patch);
children_of_function = cell (numel(ind), 1);

if (is_scalar)
  aux = cell (ndim, 1);
  for iptc = 1:npatch
    gnum = hspace.space_of_level(lev).gnum{iptc};
    [~,indices,~] = intersect (gnum, ind);

    [ind_sub{:}] = ind2sub ([hspace.space_of_level(lev).sp_patch{iptc}.ndof_dir, 1], indices); % The extra 1 makes it work in any dimension

    gnum = hspace.space_of_level(lev+1).gnum{iptc};
    for ii = 1:numel(ind_sub{1})
      for idim = 1:ndim
        aux{idim} = find (hspace.Proj{lev, iptc}{idim}(:,ind_sub{idim}(ii)));
      end
      [z{1:ndim}] = ndgrid (aux{:});
      auxI = sub2ind ([hspace.space_of_level(lev+1).sp_patch{iptc}.ndof_dir, 1], z{:});
      children = union (children, gnum(auxI(:)));
      children_of_function{ii} = gnum(auxI(:)');
    end
  end
else
  aux = cell (ndim, 1);
  for iptc = 1:npatch
    gnum = hspace.space_of_level(lev).gnum{iptc};
    [~,indices,~] = intersect (gnum, ind);
    gnum = hspace.space_of_level(lev+1).gnum{iptc};
    
    cumsum_ndof_coarse = hspace.space_of_level(lev).sp_patch{iptc}.cumsum_ndof;  
    cumsum_ndof_fine = hspace.space_of_level(lev+1).sp_patch{iptc}.cumsum_ndof;
    for icomp = 1:hspace.space_of_level(lev).sp_patch{iptc}.ncomp_param
      ind_comp = indices(indices>cumsum_ndof_coarse(icomp) & indices<=cumsum_ndof_coarse(icomp+1)) - cumsum_ndof_coarse(icomp);
      [ind_sub{:}] = ind2sub ([hspace.space_of_level(lev).sp_patch{iptc}.ndof_dir(icomp,:), 1], ind_comp); % The extra 1 makes it work in any dimension

      for ii = 1:numel(ind_sub{1})
        for idim = 1:ndim
          aux{idim} = find (hspace.Proj{lev, iptc}{icomp, idim}(:,ind_sub{idim}(ii)));
%           aux{idim} = find (hspace.Proj{lev, icomp, idim}(:,ind_sub{idim}(ii)));
        end
        [z{1:ndim}] = ndgrid (aux{:});
        auxI = sub2ind ([hspace.space_of_level(lev+1).sp_patch{iptc}.ndof_dir(icomp,:), 1], z{:});
        children = union (children, gnum(auxI(:)+cumsum_ndof_fine(icomp)));
        children_of_function{ii} = gnum(auxI(:)'+cumsum_ndof_fine(icomp));
      end
    end
  end
end

if (nargout == 2)
  flag = all (ismember (ind, hspace.active{lev}));
  if (~flag)
    flag = 2 *all (ismember (ind, union (hspace.active{lev}, hspace.deactivated{lev})));
  end
end

end
