/**

    @file      lss_helpers.hpp
    @brief     Helpers for CUDA solvers
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_HELPERS_HPP_)
#define _LSS_HELPERS_HPP_

#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cusolverSp.h>

#include "lss_macros.hpp"

namespace lss_helpers
{

/**

    @class   real_sparse_solver_cuda_helpers
    @brief   Helper for sparse CUDA solver
    @details ~

**/
class real_sparse_solver_cuda_helpers
{
  private:
    cusolverSpHandle_t solver_handle_;
    cusparseMatDescr_t mat_desc_;

  public:
    explicit real_sparse_solver_cuda_helpers() : solver_handle_{NULL}, mat_desc_{NULL}
    {
    }

    ~real_sparse_solver_cuda_helpers()
    {
        cusolverSpDestroy(solver_handle_);
    }

    inline void initialize()
    {
        LSS_ASSERT((cusolverSpCreate(&solver_handle_) == CUSOLVER_STATUS_SUCCESS),
                   " Failed to initialize sparse solver handler");
        LSS_ASSERT(cusparseCreateMatDescr(&mat_desc_) == CUSPARSE_STATUS_SUCCESS,
                   " Failed to initialize matrix descriptor");
        LSS_ASSERT(cusparseSetMatIndexBase(mat_desc_, CUSPARSE_INDEX_BASE_ZERO) == CUSOLVER_STATUS_SUCCESS,
                   " Failure while setting indexing style");
        LSS_ASSERT(cusparseSetMatType(mat_desc_, CUSPARSE_MATRIX_TYPE_GENERAL) == CUSOLVER_STATUS_SUCCESS,
                   " Failure while setting matrix type");
    }
    inline cusolverSpHandle_t const &get_solver_handle()
    {
        return solver_handle_;
    }
    inline cusparseMatDescr_t const &get_matrix_descriptor()
    {
        return mat_desc_;
    }
};

/**

    @class   real_dense_solver_cuda_helpers
    @brief   Helper for dense CUDA solver
    @details ~

**/
class real_dense_solver_cuda_helpers
{
  private:
    cusolverDnHandle_t solver_handle_;
    cublasHandle_t cublas_handle_;

  public:
    explicit real_dense_solver_cuda_helpers() : solver_handle_{NULL}, cublas_handle_{NULL}
    {
    }

    ~real_dense_solver_cuda_helpers()
    {
        cusolverDnDestroy(solver_handle_);
        cublasDestroy_v2(cublas_handle_);
    }

    inline void initialize()
    {
        LSS_ASSERT((cusolverDnCreate(&solver_handle_) == CUSOLVER_STATUS_SUCCESS),
                   " Failed to initialize dense solver handler");
        LSS_ASSERT((cublasCreate_v2(&cublas_handle_) == CUBLAS_STATUS_SUCCESS), " Failed to initialize cublas handler");
    }

    inline cusolverDnHandle_t const &get_dense_solver_handle()
    {
        return solver_handle_;
    }
    inline cublasHandle_t const &get_cublas_handle()
    {
        return cublas_handle_;
    }
};

} // namespace lss_helpers

#endif ///_LSS_HELPERS_HPP_
