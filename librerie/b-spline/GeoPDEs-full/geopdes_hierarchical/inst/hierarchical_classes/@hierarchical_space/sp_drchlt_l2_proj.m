% SP_DRCHLT_L2_PROJ: assign the degrees of freedom of Dirichlet boundaries through an L2 projection.
%
%   [u, dofs] = sp_drchlt_l2_proj (hspace, hmsh, h, sides)
%
% INPUT:
%
%  hspace: object representing the hierarchical space of trial functions (see hierarchical_space)
%  hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%  h:      function handle to compute the Dirichlet condition
%  sides:  boundary sides on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: global numbering of the corresponding basis functions
%
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
% Copyright (C) 2015 Eduardo M. Garau
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

function [u, dofs] = sp_drchlt_l2_proj (hspace, hmsh, h, drchlt_sides)

rhs  = zeros (hspace.ndof, 1);

if (hmsh.ndim == 1)
  dofs = []; u = zeros (numel(drchlt_sides), 1);
  for ii = 1:numel(drchlt_sides)
    iside = drchlt_sides(ii);
    dofs = [dofs, hspace.boundary(iside).dofs];
    if (iside == 1)
      u(ii) = h(hmsh.mesh_of_level(1).breaks{1}(1), iside);
    else
      u(ii) = h(hmsh.mesh_of_level(1).breaks{1}(end), iside); 
    end
  end
  u = u(:);
  return
end


dofs = [];
nent = 0;
for iside = drchlt_sides
  nsh_max = hspace.boundary(iside).space_of_level(end).nsh_max;
  nent = nent + hmsh.boundary(iside).nel * nsh_max^2;
  dofs = union (dofs, hspace.boundary(iside).dofs);
end
rows = zeros (nent, 1);
cols = zeros (nent, 1);
vals = zeros (nent, 1);
    
ncounter = 0;
for iside = drchlt_sides
% Restrict the function handle to the specified side, in any dimension, hside = @(x,y) h(x,y,iside)
  hside = @(varargin) h(varargin{:},iside);
  f_one = @(varargin) ones (size(varargin{1}));
  [rs, cs, vs] = ...
     op_u_v_hier (hspace.boundary(iside), hspace.boundary(iside), hmsh.boundary(iside), f_one);
  bnd_dofs = hspace.boundary(iside).dofs;

  rows(ncounter+(1:numel(rs))) = bnd_dofs(rs);
  cols(ncounter+(1:numel(rs))) = bnd_dofs(cs);
  vals(ncounter+(1:numel(rs))) = vs;
  ncounter = ncounter + numel (rs);

  rhs(bnd_dofs) = rhs(bnd_dofs) + op_f_v_hier (hspace.boundary(iside), hmsh.boundary(iside), hside);
end

M = sparse (rows(1:ncounter), cols(1:ncounter), vals(1:ncounter));
u = M(dofs, dofs) \ rhs(dofs, 1);

end
