/**

    @file      lss_dense_solvers_cuda.hpp
    @brief     Represents CUDA dense solver
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_DENSE_SOLVERS_CUDA_HPP_)
#define _LSS_DENSE_SOLVERS_CUDA_HPP_

#include "../common/lss_enumerations.hpp"
#include "../common/lss_utility.hpp"
#include "../containers/lss_flat_matrix.hpp"
#include <thrust/host_vector.h>
#include <type_traits>

namespace lss_dense_solvers
{

using lss_containers::flat_matrix;
using lss_enumerations::factorization_enum;
using lss_enumerations::flat_matrix_sort_enum;
using lss_utility::container_t;
using lss_utility::sptr_t;
using lss_utility::uptr_t;

class real_dense_solver_cuda
{
  private:
    int matrix_rows_;
    int matrix_cols_;
    uptr_t<flat_matrix> matrix_data_ptr_;

    thrust::host_vector<double> h_matrix_values_;
    thrust::host_vector<double> h_rhs_values_;

    void populate();

    explicit real_dense_solver_cuda();

  public:
    explicit real_dense_solver_cuda(int matrix_rows, int matrix_columns);

    virtual ~real_dense_solver_cuda();

    real_dense_solver_cuda(real_dense_solver_cuda const &) = delete;
    real_dense_solver_cuda &operator=(real_dense_solver_cuda const &) = delete;
    real_dense_solver_cuda(real_dense_solver_cuda &&) = delete;
    real_dense_solver_cuda &operator=(real_dense_solver_cuda &&) = delete;

    void initialize(int matrix_rows, int matrix_columns);

    void set_rhs(std::vector<double> const &rhs);

    void set_rhs_value(std::size_t idx, double value);

    void set_flat_dense_matrix(flat_matrix matrix);

    void solve(std::vector<double> &container, factorization_enum factorization = factorization_enum::QRMethod);

    container_t const solve(factorization_enum factorization = factorization_enum::QRMethod);
};

typedef sptr_t<real_dense_solver_cuda> real_dense_solver_cuda_ptr;

} // namespace lss_dense_solvers

#endif ///_LSS_DENSE_SOLVERS_CUDA_HPP_
