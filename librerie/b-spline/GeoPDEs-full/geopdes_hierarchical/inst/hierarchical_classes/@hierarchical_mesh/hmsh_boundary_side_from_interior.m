% HMSH_BOUNDARY_SIDE_FROM_INTERIOR: create a mesh object with quadrature points on one boundary side, 
%    but taking into account the information from the interior. This is necessary when imposing 
%    Dirichlet conditions in a weak form, or to compute normal derivatives, for instance.
%
%     msh_side = hmsh_boundary_side_from_interior (msh, iside);
%
% INPUTS:
%     
%    hmsh:   mesh object (see hierarchical_mesh)
%    iside: number of the boundary side to compute, from 1 to 2*hmsh.ndim
%
% OUTPUT:
%
%     msh_side: mesh object, with quadrature points only on the chosen side (see hierarchical_mesh)
%
% Copyright (C) 2017 Rafael Vazquez
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

function hmsh_sfi = hmsh_boundary_side_from_interior (hmsh, iside)

  hmsh_sfi = hmsh;
  hmsh_sfi.nlevels = hmsh.boundary(iside).nlevels;
  hmsh_sfi.active = hmsh.boundary(iside).active;
  hmsh_sfi.deactivated = hmsh.boundary(iside).deactivated;
  hmsh_sfi.nel = hmsh.boundary(iside).nel;
  hmsh_sfi.nel_per_level = hmsh.boundary(iside).nel_per_level;
  
  for ilev = 1:hmsh.boundary(iside).nlevels
    hmsh_sfi.mesh_of_level(ilev) = ...
      msh_boundary_side_from_interior (hmsh.mesh_of_level(ilev), iside);
    hmsh.msh_lev{ilev} = [];
  end
  hmsh_sfi.mesh_of_level(hmsh_sfi.nlevels+1:end) = [];
  hmsh_sfi.msh_lev(hmsh_sfi.nlevels+1:end) = [];
  hmsh_sfi.boundary = [];
end
