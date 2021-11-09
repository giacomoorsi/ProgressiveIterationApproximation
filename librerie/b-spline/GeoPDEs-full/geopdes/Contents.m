% GeoPDEs: An Octave/Matlab research tool for IsoGeometric Analysis of PDEs
%  Version 3.2.0  17-Dec-2019
%
% To run the first examples, use the following scripts
%   geopdes_base_examples
%   geopdes_elasticity_examples
%   geopdes_fluid_examples
%   geopdes_maxwell_examples
%
% The package contains the following functions. A detailed help for each function is available typing "help <function_name>"
%
% Solve: solvers for some simple problems, that you can modify to solve your own problem
%   solve_laplace
%   solve_laplace_iso
%   solve_laplace_eig
%   solve_laplace_collocation
%   solve_bilaplace_2d_iso
%   solve_bilaplace_gradgrad_2d_iso
%   solve_adv_diff_2d
%   solve_linear_elasticity
%   solve_maxwell_eig
%   solve_maxwell_eig_mixed1
%   solve_maxwell_eig_mixed2_2d
%   solve_maxwell_src
%   solve_navier_stokes
%   solve_stokes
%   solve_kirchhoff_love_shell
%
% Geometry: read geometry files, that can be created with the NURBS package 
%   geo_nurbs
%   geo_load
%   geo_read_nurbs
%   geo_convert_format
%   geo_deform
% 
% Mesh: functions to generate and compute quadrature points in the NURBS geometry
%   msh_gauss_nodes
%   msh_set_quad_nodes
%   msh_set_breaks
%   msh_cartesian
%   msh_restrict_to_cells
%   @msh_cartesian/msh_eval_boundary_side
%   @msh_cartesian/msh_evaluate_col
%   @msh_cartesian/msh_evaluate_element_list
%   @msh_cartesian/msh_precompute
%   @msh_cartesian/msh_boundary_side_from_interior
%
% Space: functions to generate and compute discrete spline and NURBS spaces
%   sp_bspline
%   sp_nurbs
%   sp_h2_error
%   sp_h1_error
%   sp_l2_error
%   sp_hcurl_error
%   sp_eval_msh
%   sp_to_vtk_b64
%   sp_to_vtk_raw
%   sp_bspline_fluid
%   sp_bspline_1d_param
%   sp_grad_preserving_transform
%   sp_integral_preserving_transform
%   sp_vector_grad_preserving_transform
%   sp_vector_curl_preserving_transform
%   sp_vector_div_preserving_transform 
%   sp_scalar
%   @sp_scalar/sp_eval_boundary_side
%   @sp_scalar/sp_evaluate_col
%   @sp_scalar/sp_evaluate_col_param
%   @sp_scalar/sp_evaluate_element_list
%   @sp_scalar/sp_evaluate_element_list_param
%   @sp_scalar/sp_h2_error
%   @sp_scalar/sp_h1_error
%   @sp_scalar/sp_l2_error
%   @sp_scalar/sp_precompute
%   @sp_scalar/sp_precompute_param
%   @sp_scalar/sp_eval
%   @sp_scalar/sp_to_vtk
%   @sp_scalar/sp_drchlt_l2_proj
%   @sp_scalar/sp_weak_drchlt_bc_laplace
%   @sp_scalar/sp_get_basis_functions
%   @sp_scalar/sp_get_cells
%   @sp_scalar/sp_get_neighbors
%   sp_vector
%   @sp_vector/sp_eval_boundary_side
%   @sp_vector/sp_eval
%   @sp_vector/sp_evaluate_col
%   @sp_vector/sp_evaluate_col_param
%   @sp_vector/sp_evaluate_element_list
%   @sp_vector/sp_evaluate_element_list_param
%   @sp_vector/sp_h2_error
%   @sp_vector/sp_h1_error
%   @sp_vector/sp_l2_error
%   @sp_vector/sp_hcurl_error
%   @sp_vector/sp_precompute
%   @sp_vector/sp_precompute_param
%   @sp_vector/sp_to_vtk
%   @sp_vector/sp_drchlt_l2_proj
%   @sp_vector/sp_weak_drchlt_bc_stokes
%   @sp_vector/sp_get_basis_functions
%   @sp_vector/sp_get_cells
%   @sp_vector/sp_get_neighbors
%
% Operators: functions to assemble different matrices and vectors
%   op_f_v
%   op_u_v
%   op_gradu_gradv
%   op_laplaceu_laplacev
%   op_vel_dot_gradu_v
%   op_mat_stab_SUPG
%   op_rhs_stab_SUPG
%   op_curlu_curlv_2d
%   op_curlu_curlv_3d
%   op_curlv_p
%   op_div_v_q
%   op_divu_divv
%   op_eu_ev
%   op_fdotn_vdotn
%   op_fdotn_v
%   op_f_vxn_2d
%   op_f_vxn_3d
%   op_uxn_vxn_2d
%   op_uxn_vxn_3d
%   op_gradgradu_gradgradv
%   op_gradu_v_otimes_n
%   op_gradv_n_f
%   op_gradv_n_u
%   op_gradu_n_gradv_n
%   op_pn_v
%   op_su_ev
%   op_udotn_vdotn
%   op_ugradu_v
%   op_ugradu_jac_v
%   op_u_otimes_n_v_otimes_n
%   op_v_gradp
%   op_KL_shells
%   op_KL_bending_stress
%   op_KL_membrane_stress
%
% Multipatch: functions of different kinds to solve problems in domains defined by multiple NURBS patches
%   mp_geo_load
%   mp_geo_read_nurbs
%   mp_interface
%   mp_interface_vector
%   mp_interface_hcurl
%   mp_interface_hdiv
%   mp_solve_laplace
%   mp_solve_linear_elasticity
%   mp_solve_maxwell_eig
%   mp_solve_maxwell_eig_mixed1
%   mp_solve_stokes
%   mp_solve_stokes_div_conforming
%   msh_multipatch
%   @msh_multipatch/msh_evaluate_element_list
%   sp_multipatch
%   @sp_multipatch/sp_l2_error
%   @sp_multipatch/sp_h1_error
%   @sp_multipatch/sp_hcurl_error
%   @sp_multipatch/sp_to_vtk
%   @sp_multipatch/sp_evaluate_element_list
%   @sp_multipatch/sp_get_cells
%   @sp_multipatch/sp_get_basis_functions
%   @sp_multipatch/sp_get_neighbors
%   @sp_multipatch/sp_drchlt_l2_proj
%   @sp_multipatch/sp_drchlt_l2_proj_udotn
%   @sp_multipatch/sp_weak_drchlt_bc_laplace
%   @sp_multipatch/sp_weak_drchlt_bc_stokes
%
% Utilities: other functions in the package
%   collocation_csp
%   knt_derham
%   msh_to_vtk
%   newtons_method
%
