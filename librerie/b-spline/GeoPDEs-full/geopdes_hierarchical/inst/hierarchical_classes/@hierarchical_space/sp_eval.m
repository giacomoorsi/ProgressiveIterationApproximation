% SP_EVAL: Compute the value or the derivatives of a hierarchical spline function, given by its degrees of freedom, at a given set of points.
%
%   [eu, F] = sp_eval (u, hspace, geometry, pts, [option]);
%   [eu, F] = sp_eval (u, hspace, geometry, npts, [option]);
%
% INPUT:
%     
%     u:         vector of dof weights
%     hspace:    object defining the discrete space (see hierarchical_space)
%     geometry:  geometry structure (see geo_load)
%     pts:       cell array with coordinates of points along each parametric direction
%     npts:      number of points along each parametric direction
%     option:    accepted options are 'value' (default), 'gradient', 'laplacian'
%
% OUTPUT:
%
%     eu: the function evaluated at the given points 
%     F:  grid points in the physical domain, that is, the mapped points
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

function [eu, F] = sp_eval (u, hspace, geometry, npts, varargin)

  sp_lev = hspace.space_of_level(hspace.nlevels);
  C = hspace_subdivision_matrix (hspace);
  u_lev =  C{hspace.nlevels} * u;
  [eu, F] = sp_eval (u_lev, sp_lev, geometry, npts, varargin{:});
  
end