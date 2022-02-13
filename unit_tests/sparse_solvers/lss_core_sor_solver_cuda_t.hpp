/**

    @file      lss_core_sor_solver_cuda_t.hpp
    @brief     Unit test for CUDA SOR solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_CORE_SOR_SOLVER_CUDA_T)
#define _LSS_CORE_SOR_SOLVER_CUDA_T

extern void device_sor_solver_basic_test();

extern void test_device_sor_solver_basic();

extern void device_sor_solver_ode_test();

extern void test_device_sor_solver_ode();

#endif ///_LSS_CORE_SOR_SOLVER_CUDA_T
