% MATRIX_BASIS_CHANGE__: compute the subdivision matrix between two consecutive levels.
%        This method is intended to remain private.
%
% function C = matrix_basis_change__ (hspace, lev, [ind_coarse], [ind_fine])
%
% Compute the new matrices to represent functions of level "lev-1"
% as linear combinations of splines (active and inactive) of level "lev"
%
% INPUT:  
%
%   hspace: an object of the class hierarchical_space_mp
%   lev:    the level for which we compute the matrix
%   ind_coarse: column indices for which to compute the output matrix
%   ind_fine:   row indices for which to compute the output matrix
%
% OUTPUT:
%
%   C:    matrix to change basis from level lev-1 to level lev
%
% If ind_fine is not given, the output matrix has size ndof(lev) x ndof(lev-1), 
%  independently of the presence of ind_coarse, for compatibility with previous versions
% If ind_fine is given, the output matrix has size numel(ind_fine) x numel(ind_coarse)
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2018-2019 Luca Coradello
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

function C = matrix_basis_change__ (hspace, lev, ind_coarse, ind_fine)

npatch = hspace.space_of_level(1).npatch;

sp_coarse = hspace.space_of_level(lev-1);
sp_fine = hspace.space_of_level(lev);
Proj = hspace.Proj(lev-1,:);

if (nargin == 4)
  C = sparse (numel(ind_fine), numel(ind_coarse));
  for iptc = 1:npatch
    spc_patch = sp_coarse.sp_patch{iptc};
    spf_patch = sp_fine.sp_patch{iptc};
    
    [~,local_indices_coarse, ind_c] = intersect (sp_coarse.gnum{iptc}, ind_coarse);
    [~,local_indices_fine, ind_f] = intersect (sp_fine.gnum{iptc}, ind_fine);
    Cpatch = subdivision_matrix_two_levels__ (spc_patch, spf_patch, Proj{iptc}, local_indices_coarse, local_indices_fine);
 
    C(ind_f, ind_c) = Cpatch;
  end
elseif (nargin == 3)
  C = sparse (sp_fine.ndof, sp_coarse.ndof);
  for iptc = 1:npatch
    spc_patch = sp_coarse.sp_patch{iptc};
    spf_patch = sp_fine.sp_patch{iptc};
    [~,local_indices,~] = intersect (sp_coarse.gnum{iptc}, ind_coarse);
    Cpatch = subdivision_matrix_two_levels__ (spc_patch, spf_patch, Proj{iptc}, local_indices);
    C(sp_fine.gnum{iptc},sp_coarse.gnum{iptc}) = Cpatch;
  end

else
  C = sparse (sp_fine.ndof, sp_coarse.ndof);
  for iptc = 1:npatch
    spc_patch = sp_coarse.sp_patch{iptc};
    spf_patch = sp_fine.sp_patch{iptc};
    Cpatch = subdivision_matrix_two_levels__ (spc_patch, spf_patch, Proj{iptc});
    C(sp_fine.gnum{iptc},sp_coarse.gnum{iptc}) = Cpatch;
  end
end

if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  if (nargin == 4)
    [~, indices] = ismember(indices,ind_fine);
  end
  C(indices,:) = 0;
end

end
