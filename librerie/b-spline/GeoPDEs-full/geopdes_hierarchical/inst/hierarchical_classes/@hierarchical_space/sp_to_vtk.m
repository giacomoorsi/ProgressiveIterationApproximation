% SP_TO_VTK: Export to VTK format for plotting.
%
%  sp_to_vtk (u, hspace, geometry, npts, filename, fieldname, [options], [lambda_lame, mu_lame])
%  sp_to_vtk (u, hspace, geometry, pts, filename, fieldname, [options], [lambda_lame, mu_lame])
%
% INPUT:
%     
%     u:           vector of dof weights
%     hspace:      object representing the space of discrete functions (see hierarchical_space)
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
% Copyright (C) 2015, 2016 Rafael Vazquez
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
  C = hspace_subdivision_matrix (hspace);
  u_lev =  C{hspace.nlevels} * u;
  sp_to_vtk (u_lev, sp_lev, geometry, npts, filename, fieldname, varargin{:});

% nopts = numel (varargin);
% 
% % For vector-valued spaces, the value of catdir is then corrected by adding one
% value = false; gradient = false; laplacian = false;
% hessian = false; curl = false; divergence = false;
% for iopt = 1:nopts
%     switch (lower (varargin{iopt}))
%       case 'value'
%         value = true;
%         catdir(iopt) = 2;
%       case 'gradient'
%         gradient = true;
%         catdir(iopt) = 3;
%       case 'laplacian' % Only for scalars, at least for now
%         laplacian = true;
%         catdir(iopt) = 2;
%       case 'hessian'
%         hessian = true;
%         catdir(iopt) = 4;
%       case 'curl' % Only for vectors
%         curl = true;
%         if (hspace.space_of_level(1).ncomp_param == 2)
%           catdir(iopt) = 1;
%         elseif (hspace.space_of_level(1).ncomp_param == 3)
%           catdir(iopt) = 2;
%         end
%       case 'divergence' % Only for vectors
%         divergence = true;
%         catdir(iopt) = 1;
%       otherwise
%         error ('hspace_eval_hmsh: unknown option: %s', varargin{iopt})
%     end
% end
% if (hspace.ncomp ~= 1)
%     catdir = catdir + 1;
%     eval_element_list = @(SP, MSH) sp_evaluate_element_list (SP, MSH, ...
%         'value', value, 'gradient', gradient, 'hessian', hessian, 'curl', curl, 'divergence', divergence);
% else
%     eval_element_list = @(SP, MSH) sp_evaluate_element_list (SP, MSH, ...
%         'value', value, 'gradient', gradient, 'laplacian', laplacian, 'hessian', hessian);
% end
% 
% F = []; eu = cell (nopts, 1);
% 
% last_dof = cumsum (hspace.ndof_per_level);
% 
% for iLev = 1:hspace.nlevels
%     sp_lev = hspace.space_of_level(iLev);
%     
%     u_lev =  hspace.Csub{iLev} * u(1:last_dof(iLev));
%     
%     [eu, F] = sp_eval (u, sp_lev, geometry, npts, varargin{:});
% 
%     msh_to_vtk (F, eu, filename, fieldname);
% end
end
