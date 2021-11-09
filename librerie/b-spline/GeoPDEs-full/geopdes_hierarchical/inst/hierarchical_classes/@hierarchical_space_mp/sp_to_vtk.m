% SP_TO_VTK: Export to VTK format for plotting.
%
%  sp_to_vtk (u, hspace, geometry, npts, filename, fieldname, [options], [lambda_lame, mu_lame])
%  sp_to_vtk (u, hspace, geometry, pts, filename, fieldname, [options], [lambda_lame, mu_lame])
%
% INPUT:
%     
%     u:           vector of dof weights
%     hspace:      object representing the space of discrete functions (see hierarchical_space_mp)
%     geometry:    geometry structure (see geo_load)
%     npts:        number of points along each parametric direction where to evaluate
%     pts:         cell array with the coordinates along each parametric direction of the points where to evaluate
%     filename:    name of the output file. 
%     fieldnames:  how to name the saved variables in the vtk file
%     options:     cell array with the fields to plot
%                   accepted options are 'value' (default), 'gradient',
%                   and for vectors also 'curl', 'divergence', 'stress'
%     lambda_lame: function handle to the first Lame coefficient (only needed to compute 'stress')
%     mu_lame:     function handle for the second Lame coefficient (only needed to compute 'stress')
%
% OUTPUT:
%
%    none
%
%    The current version is very unefficient, as it passes the solution to
%     the tensor-product space of the finest level.
%
% Copyright (C) 2015 Rafael Vazquez
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

function sp_to_vtk (u, hspace, geometry, npts, filename, fieldname, varargin)

  sp_lev = hspace.space_of_level(hspace.nlevels);
  C = hspace_subdivision_matrix (hspace, [], 'full');
  u_lev =  C{hspace.nlevels} * u;
  sp_to_vtk (u_lev, sp_lev, geometry, npts, filename, fieldname, varargin{:});

end
