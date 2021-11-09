% OP_F_V_HIER: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i) for hierarchical splines, 
%  exploiting the multilevel structure.
%
%   rhs = op_f_v_hier (hspace, hmsh, coeff, [patches]);
%
% INPUT:
%     
%   hspace: object representing the hierarchical space of test functions (see hierarchical_space_mp)
%   hmsh:   object representing the hierarchical mesh (see hierarchical_mesh_mp)
%   coeff:  function handle to compute the source function
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   rhs: assembled right-hand side
% 
% The multilevel structure is exploited in such a way that only functions
%  of the same level of the active elements have to be computed. See also
%  op_gradu_gradv_hier for more details.
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
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

function rhs = op_f_v_hier (hspace, hmsh, f, patch_list)

  if (nargin < 4)
    patch_list = 1:hmsh.mesh_of_level(1).npatch;
  end

  rhs = zeros (hspace.ndof, 1);

  ndofs = 0;
  for ilev = 1:hmsh.nlevels
    ndofs = ndofs + hspace.ndof_per_level(ilev);
    if (hmsh.nel_per_level(ilev) > 0)
      msh_lev = msh_restrict_to_patches (hmsh.msh_lev{ilev}, patch_list);

      if (msh_lev.nel > 0)
        x = cell (hmsh.rdim, 1);
        for idim = 1:hmsh.rdim
          x{idim} = reshape (msh_lev.geo_map(idim,:,:), msh_lev.nqn, msh_lev.nel);
        end
        sp_lev = sp_evaluate_element_list (hspace.space_of_level(ilev), msh_lev, 'value', true);
        
        sp_lev = change_connectivity_localized_Csub (sp_lev, hspace, ilev);
        
        b_lev = op_f_v (sp_lev, msh_lev, f(x{:}));

        dofs = 1:ndofs;
        rhs(dofs) = rhs(dofs) + hspace.Csub{ilev}.' * b_lev;
      end
    end
  end

end
