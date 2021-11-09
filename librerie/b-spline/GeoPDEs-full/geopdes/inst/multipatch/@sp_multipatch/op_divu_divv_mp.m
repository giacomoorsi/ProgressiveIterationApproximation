% OP_DIVU_DIVV_MP: assemble the matrix A = [a(i,j)], a(i,j) = (epsilon div u_j, div v_i), in a multipatch domain.
%
%   mat = op_divu_divv_mp (spu, spv, msh, [epsilon], [patches])
%
% INPUT:
%
%   spu:     object that defines the space of trial functions (see sp_multipatch)
%   spv:     object that defines the space of test functions (see sp_multipatch)
%   msh:     object that defines the domain partition and the quadrature rule (see msh_multipatch)
%   epsilon: function handle to compute some physical coefficient. Equal to one if left empty.
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat:    assembled matrix
% 
% Copyright (C) 2015, 2016, 2017 Rafael Vazquez
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

function A = op_divu_divv_mp (spu, spv, msh, coeff, patch_list)

  if (nargin < 5)
    patch_list = 1:msh.npatch;
  end

  if ((spu.npatch ~= spv.npatch) || (spu.npatch ~= msh.npatch))
    error ('op_divu_divv_mp: the number of patches does not coincide')
  end
  
  ncounter = 0;
  for iptc = patch_list
    if (nargin < 4 || isempty (coeff))
      [rs, cs, vs] = op_divu_divv_tp (spu.sp_patch{iptc}, spv.sp_patch{iptc}, msh.msh_patch{iptc});
    else
      [rs, cs, vs] = op_divu_divv_tp (spu.sp_patch{iptc}, spv.sp_patch{iptc}, msh.msh_patch{iptc}, coeff);
    end
    rows(ncounter+(1:numel (rs))) = spv.gnum{iptc}(rs);
    cols(ncounter+(1:numel (rs))) = spu.gnum{iptc}(cs);

    if (~isempty (spv.dofs_ornt))
      vs = spv.dofs_ornt{iptc}(rs)' .* vs;
    end
    if (~isempty (spu.dofs_ornt))
      vs = vs .* spu.dofs_ornt{iptc}(cs)';
    end
    
    vals(ncounter+(1:numel (rs))) = vs;
    ncounter = ncounter + numel (rs);
  end

  A = sparse (rows, cols, vals, spv.ndof, spu.ndof);
  clear rows cols vals rs cs vs

end