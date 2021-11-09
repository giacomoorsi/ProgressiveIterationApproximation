% HSPACE_EVAL_HMSH: Compute the value or the derivatives of a hierarchical spline function, 
%  given by its degrees of freedom, at the points of the corresponding hierarchical mesh.
%
%   [eu, F] = hspace_eval_hmsh (u, hspace, hmsh, [option]);
%
% INPUT:
%     
%     u:         vector of dof weights
%     hspace:    object defining the discrete space (see hierarchical_space)
%     hmsh:      object representing the hierarchical mesh (see hierarchical_mesh)
%     option:    accepted options are 'value' (default), 'gradient', 'laplacian'
%
% OUTPUT:
%
%     eu: the function evaluated at the points in the hierarchical mesh
%     F:  grid points in the physical domain, that is, the mapped points
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

function [eu, F] = hspace_eval_hmsh (u, hspace, hmsh, options)

  output_cell = true;  
  if (nargin == 3)
    options = {'value'};
    output_cell = false;  
  elseif (~iscell (options))
    options = {options};
    output_cell = false;
  end
  nopts = numel (options);

% For vector-valued spaces, the value of catdir is then corrected by adding one
  value = false; gradient = false; laplacian = false; bilaplacian = false;
  hessian = false; curl = false; divergence = false;
  for iopt = 1:nopts
    switch (lower (options{iopt}))
      case 'value'
        value = true;
        catdir(iopt) = 2;
      case 'gradient'
        gradient = true;
        catdir(iopt) = 3;
      case 'laplacian' % Only for scalars, at least for now
        laplacian = true;
        catdir(iopt) = 2;
      case 'bilaplacian'
        bilaplacian = true;
        catdir(iopt) = 2;          
      case 'hessian'
        hessian = true;
        catdir(iopt) = 4;
      case 'curl' % Only for vectors
        curl = true;
        if (hspace.space_of_level(1).ncomp_param == 2)
          catdir(iopt) = 1;
        elseif (hspace.space_of_level(1).ncomp_param == 3)
          catdir(iopt) = 2;
        end
      case 'divergence' % Only for vectors
        divergence = true;
        catdir(iopt) = 1;
      otherwise
        error ('hspace_eval_hmsh: unknown option: %s', options{iopt})
    end
  end
  if (hspace.ncomp ~= 1)
    catdir = catdir + 1;
    eval_element_list = @(SP, MSH) sp_evaluate_element_list (SP, MSH, ...
        'value', value, 'gradient', gradient, 'hessian', hessian, 'curl', curl, 'divergence', divergence);
  else
      % for higher order derivatives
%     eval_element_list = @(SP, MSH) sp_evaluate_element_list (SP, MSH, ...
%         'value', value, 'gradient', gradient, 'laplacian', laplacian, 'hessian', hessian, 'bilaplacian', bilaplacian);
    eval_element_list = @(SP, MSH) sp_evaluate_element_list (SP, MSH, ...
        'value', value, 'gradient', gradient, 'laplacian', laplacian, 'hessian', hessian);
  end
  eval_fun = @(U, SP, MSH) sp_eval_msh (U, SP, MSH, options);

  F = []; eu = cell (nopts, 1);

  last_dof = cumsum (hspace.ndof_per_level);
  for ilev = 1:hmsh.nlevels % Active levels
    if (hmsh.nel_per_level(ilev) > 0)
      msh_level = hmsh.msh_lev{ilev};
      sp_level = eval_element_list (hspace.space_of_level(ilev), msh_level);
      
      sp_level = change_connectivity_localized_Csub (sp_level, hspace, ilev);
      u_lev = hspace.Csub{ilev}*u(1:last_dof(ilev));
      
      [eu_lev, F_lev] = eval_fun (u_lev, sp_level, msh_level);
      for iopt = 1:nopts
        eu{iopt} = cat (catdir(iopt), eu{iopt}, eu_lev{iopt});
      end
      F = cat (3, F, F_lev);
    end
  end
  
  if (~output_cell)
    eu = eu{1};
  end

end