/**

    @file      lss_dense_solvers_policy.hpp
    @brief     Represents CUDA dense solver policies
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_DENSE_SOLVERS_POLICY_HPP_)
#define _LSS_DENSE_SOLVERS_POLICY_HPP_
#pragma warning(disable : 4267)

#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <type_traits>

namespace lss_dense_solvers_policy
{

struct dense_solver_device
{
};

/**
    @struct dense_solver_qr
    @brief  Dense QR factorization
**/
struct dense_solver_qr : public dense_solver_device
{
    static void solve(cusolverDnHandle_t cusolver_handle, cublasHandle_t cublas_handle, std::size_t n,
                      const double *d_Acopy, std::size_t lda, const double *d_b, double *d_x);
};

/**
    @struct dense_solver_qr
    @brief  Dense LU factorization
**/
struct dense_solver_lu : public dense_solver_device
{
    static void solve(cusolverDnHandle_t cusolver_handle, cublasHandle_t cublas_handle, std::size_t n,
                      const double *d_Acopy, std::size_t lda, const double *d_b, double *d_x);
};

/**
    @struct dense_solver_qr
    @brief  Dense Cholesky factorization
**/
struct dense_solver_cholesky : public dense_solver_device
{

    static void solve(cusolverDnHandle_t cusolver_handle, cublasHandle_t cublas_handle, std::size_t n,
                      const double *d_Acopy, std::size_t lda, const double *d_b, double *d_x);
};

} // namespace lss_dense_solvers_policy

#endif ///_LSS_DENSE_SOLVERS_POLICY_HPP_
