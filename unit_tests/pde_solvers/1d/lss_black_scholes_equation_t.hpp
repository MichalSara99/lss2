/**

    @file      lss_black_scholes_equation_t.hpp
    @brief     Unit tests for heat solver - Balck-Scholes equation
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_BLACK_SCHOLES_EQUATION_T_HPP_)
#define _LSS_BLACK_SCHOLES_EQUATION_T_HPP_

// ///////////////////////////////////////////////////////////////////////////
//							BLACK SCHOLES PROBLEMS
// ///////////////////////////////////////////////////////////////////////////

// ===========================================================================
// ========================== IMPLICIT SOLVERS ===============================
// ===========================================================================

// Dirichlet boundaries:
extern void impl_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_euler();

extern void impl_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson();

extern void test_impl_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr();

extern void impl_black_scholes_equation_dirichlet_bc_sor_solver_device_euler();

extern void impl_black_scholes_equation_dirichlet_bc_sor_solver_device_crank_nicolson();

extern void test_impl_black_scholes_equation_dirichlet_bc_sor_solver_device();

extern void impl_black_scholes_equation_dirichlet_bc_sor_solver_host_euler();

extern void impl_black_scholes_equation_dirichlet_bc_sor_solver_host_crank_nicolson();

extern void test_impl_black_scholes_equation_dirichlet_bc_sor_solver_host();

extern void impl_black_scholes_equation_dirichlet_bc_double_sweep_solver_euler();

extern void impl_black_scholes_equation_dirichlet_bc_double_sweep_solver_crank_nicolson();

extern void test_impl_black_scholes_equation_dirichlet_bc_double_sweep_solver();

extern void impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_euler();

extern void impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson();

extern void test_impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver();

// forward starting call:
extern void impl_fwd_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_euler();

extern void impl_fwd_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson();

extern void test_impl_fwd_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr();

// Using stepping = getting the whole surface
extern void impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_euler_stepping();

extern void impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson_stepping();

extern void test_impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_stepping();

// ===========================================================================
// ========================== EXPLICIT SOLVERS ===============================
// ===========================================================================

// Dirichlet boundaries:

extern void expl_black_scholes_equation_dirichlet_bc_barakat_clark();

extern void expl_black_scholes_equation_dirichlet_bc_saulyev();

extern void expl_black_scholes_equation_dirichlet_bc_euler();

extern void test_expl_black_scholes_equation_dirichlet_bc_ade();

#endif //_LSS_BLACK_SCHOLES_EQUATION_T_HPP_
