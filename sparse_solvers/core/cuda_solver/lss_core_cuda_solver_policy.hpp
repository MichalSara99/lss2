/**

    @file      lss_core_cuda_solver_policy.hpp
    @brief     Policies for CUDA solver
    @details   ~
    @author    Michal Sara
    @date      12.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_CORE_CUDA_SOLVER_POLICY_HPP_)
#define _LSS_CORE_CUDA_SOLVER_POLICY_HPP_

#include <cusolverSp.h>

#include <type_traits>

#include "../../../common/lss_macros.hpp"

namespace lss_core_cuda_solver_policy
{

/* Base for sparse factorization on Host */

struct sparse_solver_host
{
};

/* Sparse QR factorization on Host */
struct sparse_solver_host_qr : public sparse_solver_host
{

    static void solve(cusolverSpHandle_t handle, cusparseMatDescr_t mat_desc, int const system_size,
                      int const non_zero_size, double const *h_mat_vals, int const *h_row_counts,
                      int const *h_col_indices, double const *h_rhs, double tol, int reorder, double *h_solution,
                      int *singularity);
};

/* Sparse LU factorization on Host */

struct sparse_solver_host_lu : public sparse_solver_host
{
    static void solve(cusolverSpHandle_t handle, cusparseMatDescr_t mat_desc, int const system_size,
                      int const non_zero_size, double const *h_mat_vals, int const *h_row_counts,
                      int const *h_col_indices, double const *h_rhs, double tol, int reorder, double *h_solution,
                      int *singularity);
};

/* Sparse Cholesky factorization on Host */
struct sparse_solver_host_cholesky : public sparse_solver_host
{

    static void solve(cusolverSpHandle_t handle, cusparseMatDescr_t mat_desc, int const system_size,
                      int const non_zero_size, double const *h_mat_vals, int const *h_row_counts,
                      int const *h_col_indices, double const *h_rhs, double tol, int reorder, double *h_solution,
                      int *singularity);
};

/* Base for sparse factorization on Device */

struct sparse_solver_device
{
};

/* Sparse QR factorization on Device */

struct sparse_solver_device_qr : public sparse_solver_device
{

    static void solve(cusolverSpHandle_t handle, cusparseMatDescr_t mat_desc, int const system_size,
                      int const non_zero_size, double const *d_mat_vals, int const *d_row_counts,
                      int const *d_col_indices, double const *d_rhs, double tol, int reorder, double *d_solution,
                      int *singularity);
};

/* Sparse Cholesky factorization on Device */

struct sparse_solver_device_cholesky : public sparse_solver_device
{

    static void solve(cusolverSpHandle_t handle, cusparseMatDescr_t mat_desc, int const system_size,
                      int const non_zero_size, double const *d_mat_vals, int const *d_row_counts,
                      int const *d_col_indices, double const *d_rhs, double tol, int reorder, double *d_solution,
                      int *singularity);
};

} // namespace lss_core_cuda_solver_policy

#endif ///_LSS_CORE_CUDA_SOLVER_POLICY_HPP_
