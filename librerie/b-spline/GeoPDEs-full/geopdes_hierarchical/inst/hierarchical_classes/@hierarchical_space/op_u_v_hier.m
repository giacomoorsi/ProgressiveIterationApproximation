% OP_U_V_HIER: assemble the mass matrix A = [a(i,j)], a(i,j) = (mu u_j, v_i), 
%  for hierarchical splines, exploiting the multilevel structure.
%
%   mat = op_u_v_hier (hspu, hspv, hmsh, [coeff]);
%   [rows, cols, values] = op_u_v_hier (hspu, hspv, hmsh, [coeff]);
%
% INPUT:
%
%   hspu:  object representing the hierarchical space of trial functions (see hierarchical_space)
%   hspv:  object representing the hierarchical space of test functions  (see hierarchical_space)
%   hmsh:  object representing the hierarchical mesh (see hierarchical_mesh)
%   coeff: function handle to compute the reaction coefficient
%
% OUTPUT:
%
%   mat:    assembled mass matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% The multilevel structure is exploited in such a way that only functions
%  of the same level of the active elements have to be computed. See also
%  op_gradu_gradv_hier for more details.
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

function varargout = op_u_v_hier (hspu, hspv, hmsh, coeff)

  if (nargin == 3)
    coeff = @(varargin) ones (size(varargin{1}));
  end
  
  M = spalloc (hspv.ndof, hspu.ndof, 3*hspu.ndof);
  
  ndofs_u = 0;
  ndofs_v = 0;
  for ilev = 1:hmsh.nlevels
    ndofs_u = ndofs_u + hspu.ndof_per_level(ilev);
    ndofs_v = ndofs_v + hspv.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      x = cell(hmsh.msh_lev{ilev}.rdim,1);
      for idim = 1:hmsh.rdim
        x{idim} = reshape (hmsh.msh_lev{ilev}.geo_map(idim,:,:), hmsh.mesh_of_level(ilev).nqn, hmsh.nel_per_level(ilev));
      end
      spu_lev = sp_evaluate_element_list (hspu.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);
      spv_lev = sp_evaluate_element_list (hspv.space_of_level(ilev), hmsh.msh_lev{ilev}, 'value', true);

      spu_lev = change_connectivity_localized_Csub (spu_lev, hspu, ilev);  
      spv_lev = change_connectivity_localized_Csub (spv_lev, hspv, ilev);
      
      M_lev = op_u_v (spu_lev, spv_lev, hmsh.msh_lev{ilev}, coeff (x{:}));

      dofs_u = 1:ndofs_u;
      dofs_v = 1:ndofs_v;
      M(dofs_v,dofs_u) = M(dofs_v,dofs_u) + hspv.Csub{ilev}.' * M_lev * hspu.Csub{ilev};
    end
  end
  
  if (nargout == 1)
    varargout{1} = M;
  elseif (nargout == 3)
    [rows, cols, vals] = find (M);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end
end
