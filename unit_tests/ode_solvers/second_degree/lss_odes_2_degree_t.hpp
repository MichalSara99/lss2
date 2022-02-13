/**

    @file      lss_odes_2_degree_t.hpp
    @brief     Unit tests for 2nd degree ODEs
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_ODE_2_DEGREE_T_HPP_)
#define _LSS_ODE_2_DEGREE_T_HPP_

extern void impl_simple_ode_dirichlet_bc_cuda_solver_device_qr();

extern void test_impl_simple_ode_dirichlet_bc_cuda_solver_device();

extern void impl_simple_ode_dirichlet_neumann_bc_cuda_solver_device_qr();

extern void test_impl_simple_ode_dirichlet_neumann_bc_cuda_solver_device();

extern void impl_simple_ode_dirichlet_robin_bc_cuda_solver_device_qr();

extern void test_impl_simple_ode_dirichlet_robin_bc_cuda_solver_device();

extern void impl_simple_ode_neumann_robin_bc_cuda_solver_device_qr();

extern void test_impl_simple_ode_neumann_robin_bc_cuda_solver_device();

extern void impl_simple_ode1_neumann_robin_bc_cuda_solver_device_qr();

extern void test_impl_simple_ode1_neumann_robin_bc_cuda_solver_device();

#endif ///_LSS_ODE_2_DEGREE_T_HPP_
