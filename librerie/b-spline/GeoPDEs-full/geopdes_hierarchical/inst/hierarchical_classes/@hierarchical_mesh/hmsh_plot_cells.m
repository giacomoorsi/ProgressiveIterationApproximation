% HMSH_PLOT_CELLS: plot the cells of the hierarchical mesh.
%
%   hmsh_plot_cells (hmsh, [npts, fig_number])
%
% INPUT:
%
%    hmsh:       object representing the hierarchical mesh (see hierarchical_mesh)
%    npts:       number of points to use on each edge of the cell
%    fig_number: figure number where to plot (if not given, plot in the current figure)
%
% To plot the cells of the mesh, the function plots its edges using plot3.
%  A relative big number of points is needed for curved geometries.
% The function is still unefficient, but it can be useful
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

function nfig = hmsh_plot_cells (hmsh, npts, nfig)

if (nargin == 3 && ishandle (nfig))
  set (groot, 'CurrentFigure', nfig);
elseif (nargout == 1)
  nfig = gcf;
end

if (nargin < 2 || isempty (npts))
  npts = 20;
end
hold_flag = ishold;


for ilev = 1:hmsh.nlevels
  if (hmsh.nel_per_level(ilev) > 0)
    rule = cell(hmsh.ndim,1);
    for idim = 1:hmsh.ndim
      rule{idim} = [linspace(-1+1e-12, 1-1e-12, npts); zeros(1, npts)];
    end
    aux_geometry.rdim = hmsh.rdim;
    aux_geometry.map = hmsh.mesh_of_level(ilev).map;
    aux_geometry.map_der = hmsh.mesh_of_level(ilev).map_der;
    qn = msh_set_quad_nodes (hmsh.mesh_of_level(ilev).breaks, rule);
    msh_plot = msh_cartesian (hmsh.mesh_of_level(ilev).breaks, qn, [], aux_geometry, 'boundary', false);
    msh_level = msh_evaluate_element_list (msh_plot, hmsh.active{ilev});
    clear msh_plot

    x = cell (hmsh.rdim, 1);
    for idim = 1:hmsh.rdim
      x{idim} = reshape (msh_level.geo_map(idim,:,:), [npts * ones(1, hmsh.ndim), msh_level.nel]);
    end
    for idim = hmsh.rdim+1:3
      x{idim} = zeros (size (x{1}));
    end
    
    if (hmsh.ndim == 1)
      plot3 (x{1}([1 end],:), x{2}([1 end],:), x{3}([1 end],:), 'k', 'Marker', 'x');
      hold on
      plot3 (x{1}, x{2}, x{3}, 'k-');
    elseif (hmsh.ndim == 2)
      for iel = 1:msh_level.nel
        plot3 (x{1}(1,:,iel), x{2}(1,:,iel), x{3}(1,:,iel), 'k');
        hold on
        plot3 (x{1}(end,:,iel), x{2}(end,:,iel), x{3}(end,:,iel), 'k');
        plot3 (x{1}(:,1,iel), x{2}(:,1,iel), x{3}(:,1,iel), 'k');
        plot3 (x{1}(:,end,iel), x{2}(:,end,iel), x{3}(:,end,iel), 'k');
      end
    elseif (hmsh.ndim == 3)
      for iel = 1:msh_level.nel
        plot3 (x{1}(1,1,:,iel), x{2}(1,1,:,iel), x{3}(1,1,:,iel), 'k');
        hold on
        plot3 (x{1}(end,1,:,iel), x{2}(end,1,:,iel), x{3}(end,1,:,iel), 'k');
        plot3 (x{1}(1,end,:,iel), x{2}(1,end,:,iel), x{3}(1,end,:,iel), 'k');
        plot3 (x{1}(end,end,:,iel), x{2}(end,end,:,iel), x{3}(end,end,:,iel), 'k');
        plot3 (x{1}(1,:,1,iel), x{2}(1,:,1,iel), x{3}(1,:,1,iel), 'k');
        plot3 (x{1}(end,:,1,iel), x{2}(end,:,1,iel), x{3}(end,:,1,iel), 'k');
        plot3 (x{1}(1,:,end,iel), x{2}(1,:,end,iel), x{3}(1,:,end,iel), 'k');
        plot3 (x{1}(end,:,end,iel), x{2}(end,:,end,iel), x{3}(end,:,end,iel), 'k');
        plot3 (x{1}(:,1,1,iel), x{2}(:,1,1,iel), x{3}(:,1,1,iel), 'k');
        plot3 (x{1}(:,end,1,iel), x{2}(:,end,1,iel), x{3}(:,end,1,iel), 'k');
        plot3 (x{1}(:,1,end,iel), x{2}(:,1,end,iel), x{3}(:,1,end,iel), 'k');
        plot3 (x{1}(:,end,end,iel), x{2}(:,end,end,iel), x{3}(:,end,end,iel), 'k');
      end
    end
  end
end


if (~hold_flag)
  hold off
end

end
