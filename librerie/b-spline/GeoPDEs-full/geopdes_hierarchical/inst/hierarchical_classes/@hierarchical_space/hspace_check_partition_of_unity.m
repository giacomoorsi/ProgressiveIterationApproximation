% HSPACE_CHECK_PARTITION_OF_UNITY: checks the partition of unity at the
%  points in hmsh, using the coefficients in hspace.coeff_pou.
%
%    value = hspace_check_partition_of_unity (hspace, hmsh);
%    value = hspace_check_partition_of_unity (hspace, hmsh, tol);
%
% INPUT:
%
%    hspace: object representing the hierarchical space (see hierarchical_space)
%    hmsh:   object representing the hierarchical mesh (see hierarchical_mesh)
%    tol:    a given tolerance. The default value is 1e-12.
%
% OUTPUT:
%
%    value:  true if the partition of unity is satisfied, false otherwise.
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
%

function value = hspace_check_partition_of_unity (hspace, hmsh, tol)

if (nargin < 3)
  tol = 1e-10;
end

Z = hspace_eval_hmsh (hspace.coeff_pou, hspace, hmsh);

aux = max(abs(Z(:)-1));
value = (aux < tol);

if (value == 0)
    fprintf('Warning: max(abs(Z(:)-1)) = %g\n', aux);
end

end