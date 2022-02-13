/**

    @file      lss_pure_heat_equation_t.hpp
    @brief     Unit tests for heat solver - pure heat equation
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_PURE_HEAT_EQUATION_T_HPP_)
#define _LSS_PURE_HEAT_EQUATION_T_HPP_

// ///////////////////////////////////////////////////////////////////////////
//							PURE HEAT PROBLEMS
// ///////////////////////////////////////////////////////////////////////////

// ===========================================================================
// ========================== IMPLICIT SOLVERS ===============================
// ===========================================================================

// ===========================================================================
// =========== Heat problem with homogeneous boundary conditions =============
// ===========================================================================

// Dirichlet boundaries:

extern void impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_euler();

extern void impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson();

extern void test_impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr();

extern void impl_pure_heat_equation_dirichlet_bc_sor_solver_device_euler();

extern void impl_pure_heat_equation_dirichlet_bc_sor_solver_device_crank_nicolson();

extern void test_impl_pure_heat_equation_dirichlet_bc_sor_solver_device();

extern void impl_pure_heat_equation_dirichlet_bc_sor_solver_host_euler();

extern void impl_pure_heat_equation_dirichlet_bc_sor_solver_host_crank_nicolson();

extern void test_impl_pure_heat_equation_dirichlet_bc_sor_solver_host();

extern void impl_pure_heat_equation_dirichlet_bc_cuda_solver_host_qr_euler();

extern void impl_pure_heat_equation_dirichlet_bc_cuda_solver_host_qr_crank_nicolson();

extern void test_impl_pure_heat_equation_dirichlet_bc_cuda_solver_host_qr();

extern void impl_pure_heat_equation_dirichlet_bc_double_sweep_solver_euler();

extern void impl_pure_heat_equation_dirichlet_bc_double_sweep_solver_crank_nicolson();

extern void test_impl_pure_heat_equation_dirichlet_bc_double_sweep_solver();

extern void impl_pure_heat_equation_dirichlet_bc_thomas_lu_solver_euler();

extern void impl_pure_heat_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson();

extern void test_impl_pure_heat_equation_dirichlet_bc_thomas_lu_solver();

// Neuman-Dirichlet Boundaries:

extern void impl_pure_heat_equation_neumann_bc_cuda_solver_device_qr_euler();

extern void impl_pure_heat_equation_neumann_bc_cuda_solver_device_qr_crank_nicolson();

extern void test_impl_pure_heat_equation_neumann_bc_cuda_solver_device_qr();

// Neumann-Neumann Boundaries:
extern void impl_pure_heat_equation_neumann_bc_thomas_lu_solver_euler();

extern void impl_pure_heat_equation_neumann_bc_thomas_lu_solver_crank_nicolson();

extern void test_impl_pure_heat_equation_neumann_bc_thomas_lu_solver();

extern void impl_pure_heat_equation_neumann_bc_double_sweep_solver_euler();

extern void impl_pure_heat_equation_neumann_bc_double_sweep_solver_crank_nicolson();

extern void test_impl_pure_heat_equation_neumann_bc_double_sweep_solver();

// get the whole surface with stepping:
extern void impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_euler_stepping();

extern void impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson_stepping();

extern void test_impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_stepping();

// ===========================================================================
// ====== Heat problem with homogeneous boundary conditions & source =========
// ===========================================================================

extern void impl_pure_heat_equation_source_dirichlet_bc_cuda_solver_device_qr_euler();

extern void impl_pure_heat_equation_source_dirichlet_bc_cuda_solver_device_qr_clark_nicolson();

extern void test_impl_pure_heat_equation_source_dirichlet_bc_cuda_solver_device();

extern void impl_pure_heat_equation_source_dirichlet_bc_sor_solver_device_euler();

extern void impl_pure_heat_equation_source_dirichlet_bc_sor_solver_device_clark_nicolson();

extern void test_impl_pure_heat_equation_source_dirichlet_bc_sor_solver_device();

// ===========================================================================
// ========================== EXPLICIT SOLVERS ===============================
// ===========================================================================

// ===========================================================================
// =========== Heat problem with homogeneous boundary conditions =============
// ===========================================================================

// Dirichlet boundaries:

extern void expl_pure_heat_equation_dirichlet_bc_barakat_clark();

extern void expl_pure_heat_equation_dirichlet_bc_saulyev();

extern void expl_pure_heat_equation_dirichlet_bc_euler();

extern void test_expl_pure_heat_equation_dirichlet_bc_ade();

extern void empl_pure_heat_equation_neumann_dirichlet_bc_euler();

extern void expl_pure_heat_equation_neumann_neumann_bc_euler();

extern void test_expl_pure_heat_equation_neumann_bc_euler();

extern void expl_pure_heat_equation_dirichlet_bc_euler_device();

extern void test_expl_pure_heat_equation_dirichlet_bc_device();

#endif //_LSS_general_heat_equation_T_HPP_
