/**

    @file      lss_sor_solver_cuda_t.hpp
    @brief     Unit tests for CUDA SOR Solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_SOR_SOLVER_CUDA_T_HPP_)
#define _LSS_SOR_SOLVER_CUDA_T_HPP_

extern void device_sor_bvp_dirichlet_bc_test();

extern void test_device_sor_bvp_dirichlet_bc();

extern void device_sor_bvp_robin_bc_test();

extern void test_device_sor_bvp_robin_bc();

extern void device_sor_bvp_mix_bc_test();

extern void test_device_sor_bvp_mix_bc();

#endif ///_LSS_SOR_SOLVER_CUDA_T_HPP_
