/**

    @file      lss_pure_wave_equation_t.hpp
    @brief     Unit tests for wave solver - pure wave equation
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_PURE_WAVE_EQUATION_T_HPP_)
#define _LSS_PURE_WAVE_EQUATION_T_HPP_

// ///////////////////////////////////////////////////////////////////////////
//							PURE WAVE PROBLEMS
// ///////////////////////////////////////////////////////////////////////////

// ===========================================================================
// ========================== IMPLICIT SOLVERS ===============================
// ===========================================================================

// ===========================================================================
// =========== Wave problem with homogeneous boundary conditions =============
// ===========================================================================

// Dirichlet boundaries:

extern void impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_qr_detail();

extern void test_impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_qr();

extern void impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_sor_detail();

extern void test_impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_sor();

extern void impl_pure_wave_equation_dirichlet_bc_cuda_solver_host_sor_detail();

extern void test_impl_pure_wave_equation_dirichlet_bc_cuda_solver_host_sor();

extern void impl_pure_wave_equation_dirichlet_bc_solver_host_double_sweep_detail();

extern void test_impl_pure_wave_equation_dirichlet_bc_solver_host_double_sweep();

extern void impl_pure_wave_equation_dirichlet_bc_solver_host_lu_detail();

extern void test_impl_pure_wave_equation_dirichlet_bc_solver_host_lu();

extern void impl_wave_equation_dirichlet_bc_solver_host_lu_detail();

extern void test_impl_wave_equation_dirichlet_bc_solver_host_lu();

extern void impl_damped_wave_equation_dirichlet_bc_solver_host_double_sweep_detail();

extern void test_impl_damped_wave_equation_dirichlet_bc_solver_host_double_sweep();

// Neumann BC:

extern void impl_pure_wave_equation_neumann_bc_cuda_solver_device_qr_detail();

extern void test_impl_pure_wave_equation_neumann_bc_cuda_solver_device_qr();

// ===========================================================================
// ========================== EXPLICIT SOLVERS ===============================
// ===========================================================================

// ===========================================================================
// =========== Wave problem with homogeneous boundary conditions =============
// ===========================================================================

// Dirichlet boundaries:

extern void expl_pure_wave_equation_dirichlet_bc_cuda_host_solver_detail();

extern void test_expl_pure_wave_equation_dirichlet_bc_cuda_host_solver();

extern void expl_pure_wave_equation_dirichlet_bc_cuda_device_solver_detail();

extern void test_expl_pure_wave_equation_dirichlet_bc_cuda_device_solver();

#endif //_LSS_PURE_WAVE_EQUATION_T_HPP_
