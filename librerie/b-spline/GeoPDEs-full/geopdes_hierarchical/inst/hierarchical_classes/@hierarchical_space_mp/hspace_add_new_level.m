% HSPACE_ADD_NEW_LEVEL: add an empty level to the hierarchical mesh.
%
%   hspace = hspace_add_new_level (hspace, hmsh)
%
% INPUT:
%
%   hspace:  object representing the hierarchical space (see hierarchical_space_mp)
%   hmsh:    object representing the hierarchical mesh, already with the new level (see hierarchical_mesh_mp)
%
% OUTPUT:
%
%   hspace:   the object of the hierarchical space with one more level, without active functions
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

function hspace = hspace_add_new_level (hspace, hmsh)

if (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_scalar'))
  is_scalar = true;
elseif (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_vector'))
  is_scalar = false;
else
  error ('Unknown space type')
end

% Computation of a multi-patch tensor product space if a new level is activated,
%  and the 1D projectors between the previous level and the new one.
if (numel(hspace.space_of_level) == hmsh.nlevels-1)
  regularity = hspace.regularity;
  if (is_scalar)
    degree = hspace.space_of_level(hmsh.nlevels-1).sp_patch{1}.degree;
  else
    scalar_spaces = hspace.space_of_level(hmsh.nlevels-1).sp_patch{1}.scalar_spaces;
    for icomp = 1:hspace.space_of_level(hmsh.nlevels-1).sp_patch{1}.ncomp_param
      degree{icomp} = scalar_spaces{icomp}.degree;
    end
  end
  msh_level = hmsh.mesh_of_level(hmsh.nlevels);
  [new_space, Proj] = sp_refine (hspace.space_of_level(hmsh.nlevels-1), msh_level, hmsh.nsub, degree, regularity);
  hspace.space_of_level(hmsh.nlevels) = new_space; clear new_space
  hspace.Proj(hmsh.nlevels-1,:) = Proj(:);

  hspace.nlevels = hmsh.nlevels;
  hspace.active{hmsh.nlevels} = [];
  hspace.deactivated{hmsh.nlevels} = [];
  hspace.ndof_per_level(hmsh.nlevels) = 0;
else
  warning ('The number of levels of the space in input should be one less than the number of levels of the mesh')
end

end