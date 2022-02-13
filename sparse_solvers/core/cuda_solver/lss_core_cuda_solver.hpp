/**

    @file      lss_core_cuda_solver.hpp
    @brief     Cuda CORE solvers
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once

#if !defined(_LSS_CORE_CUDA_SOLVER_HPP_)
#define _LSS_CORE_CUDA_SOLVER_HPP_

#pragma warning(disable : 4267)
#pragma warning(disable : 4244)

#include <cuda_runtime.h>
#include <cusolverSp.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <type_traits>

#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_helpers.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../containers/lss_flat_matrix.hpp"

namespace lss_core_cuda_solver
{

using lss_containers::flat_matrix;
using lss_enumerations::factorization_enum;
using lss_enumerations::flat_matrix_sort_enum;
using lss_enumerations::memory_space_enum;
using lss_helpers::real_sparse_solver_cuda_helpers;
using lss_utility::container_t;
using lss_utility::copy;
using lss_utility::sptr_t;
using lss_utility::uptr_t;

template <memory_space_enum memory_space> class real_sparse_solver_cuda
{
};

/**

    @class   real_sparse_solver_cuda
    @brief   CUDAs sparse solver on HOST
    @details ~

**/
template <> class real_sparse_solver_cuda<memory_space_enum::Host>
{

  protected:
    int system_size_;
    uptr_t<flat_matrix> matrix_data_ptr_;

    thrust::host_vector<double> h_matrix_values_;
    thrust::host_vector<double> h_vector_values_; // of systemSize length
    thrust::host_vector<int> h_column_indices_;
    thrust::host_vector<int> h_row_counts_; // of systemSize + 1 length

    void build_csr();
    explicit real_sparse_solver_cuda();

  public:
    void initialize(std::size_t system_size);

    explicit real_sparse_solver_cuda(int system_size);

    virtual ~real_sparse_solver_cuda();

    real_sparse_solver_cuda(real_sparse_solver_cuda const &) = delete;
    real_sparse_solver_cuda &operator=(real_sparse_solver_cuda const &) = delete;
    real_sparse_solver_cuda(real_sparse_solver_cuda &&) = delete;
    real_sparse_solver_cuda &operator=(real_sparse_solver_cuda &&) = delete;

    std::size_t non_zero_elements() const;

    void set_rhs(container_t const &rhs);

    void set_rhs_value(std::size_t idx, double value);

    void set_flat_sparse_matrix(flat_matrix matrix);

    void solve(container_t &solution, factorization_enum factorization = factorization_enum::QRMethod);

    container_t const solve(factorization_enum factorization = factorization_enum::QRMethod);
};

using host_real_sparse_solver_cuda_ptr = sptr_t<real_sparse_solver_cuda<memory_space_enum::Host>>;

/**

    @class   real_sparse_solver_cuda
    @brief   CUDAs sparse solver on DEVICE
    @details ~

**/
template <> class real_sparse_solver_cuda<memory_space_enum::Device>
{

  protected:
    int system_size_;
    uptr_t<flat_matrix> matrix_data_ptr_;

    thrust::host_vector<double> h_matrix_values_;
    thrust::host_vector<double> h_vector_values_; // of systemSize length
    thrust::host_vector<int> h_column_indices_;
    thrust::host_vector<int> h_row_counts_; // of systemSize + 1 length

    void build_csr();
    explicit real_sparse_solver_cuda();

  public:
    explicit real_sparse_solver_cuda(int system_size);

    virtual ~real_sparse_solver_cuda();

    real_sparse_solver_cuda(real_sparse_solver_cuda const &) = delete;
    real_sparse_solver_cuda &operator=(real_sparse_solver_cuda const &) = delete;
    real_sparse_solver_cuda(real_sparse_solver_cuda &&) = delete;
    real_sparse_solver_cuda &operator=(real_sparse_solver_cuda &&) = delete;

    void initialize(std::size_t system_size);

    std::size_t non_zero_elements() const;

    void set_rhs(container_t const &rhs);

    void set_rhs_value(std::size_t idx, double value);

    void set_flat_sparse_matrix(flat_matrix matrix);

    void solve(container_t &solution, factorization_enum factorization = factorization_enum::QRMethod);

    container_t const solve(factorization_enum factorization = factorization_enum::QRMethod);
};

using device_real_sparse_solver_cuda_ptr = sptr_t<real_sparse_solver_cuda<memory_space_enum::Device>>;

} // namespace lss_core_cuda_solver

#endif ///_LSS_CORE_CUDA_SOLVER_HPP_
