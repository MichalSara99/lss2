/**

    @file      lss_thomas_lu_solver_t.hpp
    @brief     Unit tests for tridiagonal Thomas LU Solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_THOMAS_LU_SOLVER_T_HPP_)
#define _LSS_THOMAS_LU_SOLVER_T_HPP_

extern void thomas_lu_bvp_dirichlet_bc_test();

extern void test_thomas_lu_bvp_dirichlet_bc();

extern void thomas_lu_bvp_robin_bc_test();

extern void test_thomas_lu_bvp_robin_bc();

extern void thomas_lu_bvp_mix_bc_test();

extern void test_thomas_lu_bvp_mix_bc();

#endif ///_LSS_THOMAS_LU_SOLVER_T_HPP_
