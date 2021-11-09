% HSPACE_GET_PARENTS: compute the parents of a given set of functions of the same level.
%
%     [parents, flag] = hspace_get_parents (hspace, lev, ind)
%
% Get the parents of a given set of basis functions of level lev, with the
%  subdivision given by the "Proj" matrices.
% All the parent functions are stored in the same array.
%
% INPUT:
%
%     hspace: the hierarchical space (see hierarchical_space_mp)
%     lev:    level of the children functions
%     ind:    indices of the functions in the multipatch space of level lev
%
% OUTPUT:
%
%     parents:  indices of the parents, with the numbering of the tensor product space
%     flag:     a flag to tell whether all the input functions are active (1) 
%               active or deactivated (2), or if there is any passive function (0)
%
% Copyright (C) 2016 Rafael Vazquez
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

function [parents, flag] = hspace_get_parents (hspace, lev, ind)

% % This computation would not work for the truncated basis
% Cmat = matrix_basis_change__ (hspace, lev);  
% [~,parents] = find (Cmat(ind,:));
% parents = unique (parents);

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
parents = [];

npatch = numel (hspace.space_of_level(1).sp_patch);

if (is_scalar)
  aux = cell (ndim, 1);
  for iptc = 1:npatch
    gnum = hspace.space_of_level(lev).gnum{iptc};
    [~,indices,~] = intersect (gnum, ind);

    [ind_sub{:}] = ind2sub ([hspace.space_of_level(lev).sp_patch{iptc}.ndof_dir, 1], indices); % The extra 1 makes it work in any dimension

    gnum = hspace.space_of_level(lev-1).gnum{iptc};
    for ii = 1:numel(ind_sub{1})
      for idim = 1:ndim
        aux{idim} = find (hspace.Proj{lev-1, iptc}{idim}(ind_sub{idim}(ii),:));
      end
      [z{1:ndim}] = ndgrid (aux{:});
      auxI = sub2ind ([hspace.space_of_level(lev-1).sp_patch{iptc}.ndof_dir, 1], z{:});
      parents = union (parents, gnum(auxI(:)));
    end
  end
else
  aux = cell (ndim, 1);
  for iptc = 1:npatch
    gnum = hspace.space_of_level(lev).gnum{iptc};
    [~,indices,~] = intersect (gnum, ind);
    gnum = hspace.space_of_level(lev-1).gnum{iptc};

    cumsum_ndof_coarse = hspace.space_of_level(lev-1).sp_patch{iptc}.cumsum_ndof;  
    cumsum_ndof_fine = hspace.space_of_level(lev).sp_patch{iptc}.cumsum_ndof;
    for icomp = 1:hspace.space_of_level(lev).sp_patch{iptc}.ncomp_param
      ind_comp = indices(indices>cumsum_ndof_fine(icomp) & indices<=cumsum_ndof_fine(icomp+1)) - cumsum_ndof_fine(icomp);
      [ind_sub{:}] = ind2sub ([hspace.space_of_level(lev).sp_patch{iptc}.ndof_dir(icomp,:), 1], ind_comp); % The extra 1 makes it work in any dimension

      for ii = 1:numel(ind_sub{1})
        for idim = 1:ndim
          aux{idim} = find (hspace.Proj{lev-1, iptc}{icomp, idim}(:,ind_sub{idim}(ii)));
%           aux{idim} = find (hspace.Proj{lev-1, icomp, idim}(:,ind_sub{idim}(ii)));
        end
        [z{1:ndim}] = ndgrid (aux{:});
        auxI = sub2ind ([hspace.space_of_level(lev-1).sp_patch{iptc}.ndof_dir(icomp,:), 1], z{:});
        parents = union (children, gnum(auxI(:)+cumsum_ndof_coarse(icomp)));
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
