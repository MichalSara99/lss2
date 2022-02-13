/**

    @file      lss_cuda_solver_t.hpp
    @brief     Unit tests for tridiagonal CUDA solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_CUDA_SOLVER_T_HPP_)
#define _LSS_CUDA_SOLVER_T_HPP_

extern void device_bvp_dirichlet_bc_test();

extern void test_device_bvp_dirichlet_bc();

extern void device_bvp_robin_bc_test();

extern void test_device_bvp_robin_bc();

extern void device_bvp_mix_bc_test();

extern void test_device_bvp_mix_bc();

#endif ///_LSS_CUDA_SOLVER_T_HPP_
