% PLOT_NUMERICAL_AND_EXACT_SOLUTION: personal function to plot.
%
%   fig = plot_numerical_and_exact_solution (u, hspace, geometry, npts, [uex, fig])
%
% INPUT:
%     
%     u:        vector of dof weights
%     hspace:   object defining the discrete space (see hierarchical_space)
%     geometry: geometry structure (see mp_geo_load)
%     npts:     number of points along each parametric direction
%     uex:      function handle to compute the exact solution (optional)
%     fig:      handle to the figure
%
% OUTPUT:
%
%     fig:      handle to the figure
% 
% Copyright (C) 2016-2019, Eduardo Garau, Rafael Vazquez
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


function fig = plot_numerical_and_exact_solution (u, hspace, geometry, npts, uex, fig)

if (nargin == 6 && ishandle(fig))
  figure (fig)
elseif (nargout == 1)
  fig = figure;
else
  figure;
end

if (isa (hspace, 'sp_scalar'))
  sp_aux = hspace;
  is_vector = false;
elseif (isa (hspace, 'sp_vector'))
  sp_aux = hspace.scalar_spaces{1};
  is_vector = true;
elseif (isa (hspace, 'hierarchical_space'))
  if (isa (hspace.space_of_level(1), 'sp_scalar'))
    sp_aux = hspace.space_of_level(1);
    is_vector = false;
  elseif (isa (hspace.space_of_level(1), 'sp_vector'))
    sp_aux = hspace.space_of_level(1).scalar_spaces{1};
    is_vector = true;
  end
elseif (isa (hspace, 'hierarchical_space_mp'))
  if (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_scalar'))
    sp_aux = hspace.space_of_level(1).sp_patch{1};
    is_vector = false;
  elseif (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_vector'))
    sp_aux = hspace.space_of_level(1).sp_patch{1}.scalar_spaces{1};
    is_vector = true;
  end
end
ndim = numel (sp_aux.knots);
first_knot = cellfun (@(x) x(1), sp_aux.knots);
last_knot = cellfun (@(x) x(end), sp_aux.knots);

if ((nargin > 4 && ~isempty (uex)))
  if (ndim == 3)
    warning ('The plot of the exact solution for volumes is not implemented. Plotting the numerical solution')
  elseif (is_vector)
    warning ('The plot of the exact solution for vectors is not implemented. Plotting the numerical solution')
  else
    subplot (1, 2, 2)
    for idim = 1:ndim
      x{idim} = linspace (first_knot(idim), last_knot(idim), npts(idim));
    end
    for iptc = 1:numel(geometry)
      F = geometry(iptc).map (x);
      rdim = size (F, 1);
      F = reshape (F, [rdim, npts]);

      if (ndim == 1 && rdim == 1)
        plot (F, uex(F)); hold on
      elseif (ndim == 1 && rdim == 2)
        plot3 (F(1,:), F(2,:), uex(F)); hold on
      elseif (ndim == 2 && rdim == 2)
        X = squeeze (F(1,:,:)); Y = squeeze (F(2,:,:));
        surf (X, Y, uex(X, Y)); hold on
        shading interp
        xlabel('x'); ylabel('y'); zlabel('z')
      elseif (ndim == 2 && rdim == 3)
        X = squeeze (F(1,:,:)); Y = squeeze (F(2,:,:)); Z = squeeze (F(3,:,:));
        surf (X, Y, Z, uex(X, Y, Z)); hold on
        shading interp
        xlabel('x'); ylabel('y'); zlabel('z')
      else
        error('The plot of the exact solution with ndim=%d and rdim=%d is not implemented', ndim, rdim)
      end
    end
    title ('Exact solution')
    hold off
    subplot (1, 2, 1)
  end
end

sp_plot_solution (u, hspace, geometry, npts)
shading interp
title ('Numerical solution'),

if ((nargin > 4 && ~isempty (uex)))
  xl = xlim; yl = ylim; zl = zlim;
  subplot (1, 2, 2)
  xlim(xl); ylim (yl); zlim (zl);
end

end
