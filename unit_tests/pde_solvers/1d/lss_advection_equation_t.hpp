/**

    @file      lss_advection_equation_t.hpp
    @brief     Unit tests for heat solver - advection equation
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_ADVECTION_EQUATION_T_HPP_)
#define _LSS_ADVECTION_EQUATION_T_HPP_

// ///////////////////////////////////////////////////////////////////////////
//							ADVECTION PROBLEMS
// ///////////////////////////////////////////////////////////////////////////

// ===========================================================================
// ========================== IMPLICIT SOLVERS ===============================
// ===========================================================================

// ===========================================================================
// ====== Advection Diffusion problem with homogeneous boundary conditions ===
// ===========================================================================

extern void impl_ddv_diff_equation_dirichlet_bc_cuda_solver_device_qr_euler();

extern void impl_adv_diff_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson();

extern void test_impl_adv_diff_equation_dirichlet_bc_cuda_solver_device_qr();

extern void impl_adv_diff_equation_dirichlet_bc_sor_solver_device_euler();

extern void impl_adv_diff_equation_dirichlet_bc_sor_solver_device_crank_nicolson();

extern void test_impl_adv_diff_equation_dirichlet_bc_sor_solver_device();

extern void impl_adv_diff_equation_dirichlet_bc_sor_solver_host_euler();

extern void impl_adv_diff_equation_dirichlet_bc_sor_solver_host_crank_nicolson();

extern void test_impl_adv_diff_equation_dirichlet_bc_sor_solver_host();

extern void impl_adv_diff_equation_dirichlet_bc_cuda_solver_host_qr_euler();

extern void impl_adv_diff_equation_dirichlet_bc_cuda_solver_host_qr_crank_nicolson();

extern void test_impl_adv_diff_equation_dirichlet_bc_cuda_solver_host_qr();

extern void impl_adv_diff_equation_dirichlet_bc_double_sweep_solver_euler();

extern void impl_adv_diff_equation_dirichlet_bc_double_sweep_solver_crank_nicolson();

extern void test_impl_adv_diff_equation_dirichlet_bc_double_sweep_solver();

extern void impl_adv_diff_equation_dirichlet_bc_thomas_lu_solver_euler();

extern void impl_adv_diff_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson();

extern void test_impl_adv_diff_equation_dirichlet_bc_thomas_lu_solver();

#endif //_LSS_ADVECTION_EQUATION_T_HPP_
