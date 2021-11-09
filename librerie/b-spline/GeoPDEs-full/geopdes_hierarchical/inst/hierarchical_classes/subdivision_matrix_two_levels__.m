% SUBDIVION_MATRIX_TWO_LEVELS__: compute the subdivision matrix between two consecutive levels.
%        It works for scalar-valued and vector-valued functions, single patch.
%        This method is intended to remain private.
%
% function C = subdivision_matrix_two_levels__ (sp_coarse, sp_fine, Proj, [ind_coarse, ind_fine])
%
% Compute the new matrices to represent functions of level "lev-1"
% as linear combinations of splines (active and inactive) of level "lev"
%
% INPUT:  
%
%   sp_coarse:  a space object, either sp_scalar or sp_vector
%   sp_fine:    a space object, obtained by refinement of sp_fine
%   Proj:       univariate subdivision matrices, given by sp_refine
%   ind_coarse: column indices for which to compute the output matrix
%   ind_fine:   row indices for which to compute the output matrix
%
% OUTPUT:
%
%   C:  subdivision matrix to change basis from level lev-1 to level lev
%
% If ind_fine is not given, the output matrix has size ndof(lev) x ndof(lev-1), 
%  independently of the presence of ind_coarse, for compatibility with previous versions
% If ind_fine is given, the output matrix has size numel(ind_fine) x numel(ind_coarse)
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2017-2019 Rafael Vazquez
% Copyright (C) 2018-2019 Luca Coradello
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

function varargout = subdivision_matrix_two_levels__ (sp_coarse, sp_fine, Proj, ind_coarse, ind_fine)

if (isa (sp_coarse, 'sp_scalar'))
  is_scalar = true;
elseif (isa (sp_coarse, 'sp_vector'))
  is_scalar = false;
  ncomp_param = sp_coarse.ncomp_param;
else
  error ('Unknown space type')
end
ndim = size (Proj, 2);

if (is_scalar)
  if (nargin < 4)
    C = 1;
    for idim = 1:ndim
      C = kron (Proj{1,idim}, C);
    end
    if (strcmpi (sp_coarse.space_type, 'NURBS'))
      Wlev = spdiags (sp_coarse.weights(:), 0, sp_coarse.ndof, sp_coarse.ndof);
      Wlev_fine = spdiags (1./sp_fine.weights(:), 0, sp_fine.ndof, sp_fine.ndof);
      C = Wlev_fine * C * Wlev;
    end
    
  elseif (nargin == 4)
    sub_coarse = cell (ndim, 1);
    [sub_coarse{:}] = ind2sub ([sp_coarse.ndof_dir, 1], ind_coarse);
  
    rows = zeros (prod (sp_fine.degree+1)*numel(ind_coarse), 1); cols = rows; vals = rows;
    ncounter = 0;
    for ii = 1:numel(ind_coarse)
      Caux = 1;
      for idim = 1:ndim
        Caux = kron (Proj{1,idim}(:,sub_coarse{idim}(ii)), Caux);
      end
      [ir, ~, iv] = find (Caux);
      rows(ncounter+(1:numel(ir))) = ir;
      cols(ncounter+(1:numel(ir))) = ind_coarse(ii);
      vals(ncounter+(1:numel(ir))) = iv;
      ncounter = ncounter + numel (ir);
    end
    rows = rows(1:ncounter); cols = cols(1:ncounter); vals = vals(1:ncounter);
    
    if (strcmpi (sp_coarse.space_type, 'NURBS'))
      weights_coarse = sp_coarse.weights(:);
      weights_fine = sp_fine.weights(:);
      vals = vals .* weights_coarse(cols) ./ weights_fine(rows);
    end
    
  elseif (nargin == 5)
    sub_coarse = cell (ndim, 1);
    [sub_coarse{:}] = ind2sub ([sp_coarse.ndof_dir, 1], ind_coarse);
  
    rows = zeros (prod (sp_fine.degree+1)*numel(ind_coarse), 1); cols = rows; vals = rows;
    ncounter = 0;
    for ii = 1:numel(ind_coarse)
      Caux = 1;
      for idim = 1:ndim
        Caux = kron (Proj{1,idim}(:,sub_coarse{idim}(ii)), Caux);
      end
      [ir, ~, iv] = find (Caux);
      [~,IA,IB] = intersect(ir,ind_fine);
      rows(ncounter+(1:numel(IB))) = IB;
      cols(ncounter+(1:numel(IB))) = ii;
      vals(ncounter+(1:numel(IB))) = iv(IA);
      ncounter = ncounter + numel (IB);
    end
    rows = rows(1:ncounter); cols = cols(1:ncounter); vals = vals(1:ncounter);
    
    if (strcmpi (sp_coarse.space_type, 'NURBS'))
      weights_coarse = sp_coarse.weights(:);
      weights_fine = sp_fine.weights(:);
      vals = vals .* weights_coarse(cols) ./ weights_fine(rows);
    end
  else    
    error('Wrong number of input parameters!')
  end
  
else
  if (nargin < 4)
    Caux = cell (ncomp_param, 1);
    for icomp = 1:ncomp_param
      spc_scalar = sp_coarse.scalar_spaces{icomp};
      spf_scalar = sp_fine.scalar_spaces{icomp};
      Caux{icomp} = subdivision_matrix_two_levels__ (spc_scalar, spf_scalar, Proj(icomp,:));
    end
    C = blkdiag (Caux{:});
    
  elseif (nargin == 4)
    rows = []; cols = []; vals = [];
    cumsum_ndof_coarse = sp_coarse.cumsum_ndof;  
    cumsum_ndof_fine = sp_fine.cumsum_ndof;  
    
    for icomp = 1:ncomp_param
      ind_comp = ind_coarse(ind_coarse>cumsum_ndof_coarse(icomp) & ...
                            ind_coarse<=cumsum_ndof_coarse(icomp+1)) - cumsum_ndof_coarse(icomp);

      spc_scalar = sp_coarse.scalar_spaces{icomp};
      spf_scalar = sp_fine.scalar_spaces{icomp};
      [rows_c, cols_c, vals_c] = subdivision_matrix_two_levels__ (spc_scalar, spf_scalar, Proj(icomp,:), ind_comp);
      rows = [rows; rows_c+cumsum_ndof_fine(icomp)]; 
      cols = [cols; cols_c+cumsum_ndof_coarse(icomp)]; 
      vals = [vals; vals_c];
    end
    
  elseif (nargin == 5)
    rows = []; cols = []; vals = [];
    cumsum_ndof_coarse = sp_coarse.cumsum_ndof;
    cumsum_ndof_fine = sp_fine.cumsum_ndof;
    for icomp = 1:ncomp_param
      ind_comp_coarse = ind_coarse(ind_coarse>cumsum_ndof_coarse(icomp) & ...
                               ind_coarse<=cumsum_ndof_coarse(icomp+1));  
      ind_comp_fine = ind_fine(ind_fine>cumsum_ndof_fine(icomp) & ...
                             ind_fine<=cumsum_ndof_fine(icomp+1));  
      ind_scalar_coarse = ind_comp_coarse - cumsum_ndof_coarse(icomp);
      ind_scalar_fine = ind_comp_fine - cumsum_ndof_fine(icomp);

      spc_scalar = sp_coarse.scalar_spaces{icomp};
      spf_scalar = sp_fine.scalar_spaces{icomp};
      [rows_c, cols_c, vals_c] = subdivision_matrix_two_levels__ (spc_scalar, spf_scalar, Proj(icomp,:), ind_scalar_coarse, ind_scalar_fine);
      [~,rows_c] = ismember (ind_scalar_fine(rows_c)+cumsum_ndof_fine(icomp), ind_fine);
      [~,cols_c] = ismember (ind_scalar_coarse(cols_c)+cumsum_ndof_coarse(icomp), ind_coarse);
      rows = [rows; rows_c];
      cols = [cols; cols_c];
      vals = [vals; vals_c];
    end

  else    
    error('Wrong number of input parameters!')
  end
end

if (nargout == 1)  
  if (nargin == 4)
    C = sparse (rows, cols, vals, sp_fine.ndof, sp_coarse.ndof);
  elseif (nargin == 5)
    C = sparse (rows, cols, vals, numel(ind_fine), numel(ind_coarse));
  end
  varargout{1} = C;
elseif (nargout == 3)
  varargout = {rows, cols, vals};
end

end
