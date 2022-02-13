#include "lss_dense_solvers_cuda.hpp"
#include "lss_dense_solvers_policy.hpp"

#include "../common/lss_helpers.hpp"
#include "../common/lss_macros.hpp"

#include <algorithm>
#include <cusolverDn.h>
#include <thrust/device_vector.h>

namespace lss_dense_solvers
{

using lss_dense_solvers_policy::dense_solver_cholesky;
using lss_dense_solvers_policy::dense_solver_lu;
using lss_dense_solvers_policy::dense_solver_qr;
using lss_helpers::real_dense_solver_cuda_helpers;

real_dense_solver_cuda::real_dense_solver_cuda()
{
}

real_dense_solver_cuda::real_dense_solver_cuda(int matrix_rows, int matrix_columns)
    : matrix_cols_{matrix_columns}, matrix_rows_{matrix_rows}
{
}

real_dense_solver_cuda ::~real_dense_solver_cuda()
{
}

void real_dense_solver_cuda::set_rhs(std::vector<double> const &rhs)
{
    LSS_ASSERT(rhs.size() == matrix_rows_, " Right-hand side vector of the system has incorrect size.");
    thrust::copy(rhs.begin(), rhs.end(), h_rhs_values_.begin());
}

void real_dense_solver_cuda::set_rhs_value(std::size_t idx, double value)
{
    LSS_ASSERT((idx < matrix_rows_), "idx is outside range");
    h_rhs_values_[idx] = value;
}

void real_dense_solver_cuda::set_flat_dense_matrix(flat_matrix matrix)
{
    LSS_ASSERT((matrix.rows() == matrix_rows_) && (matrix.columns() == matrix_cols_),
               " flat_matrix has incorrect number of rows or columns");
    matrix_data_ptr_ = std::make_unique<flat_matrix>(std::move(matrix));
}

void real_dense_solver_cuda::populate()
{
    LSS_VERIFY(matrix_data_ptr_, "flat_matrix has not been provided.");
    // CUDA Dense solver is column-major:
    matrix_data_ptr_->sort(flat_matrix_sort_enum::ColumnMajor);

    for (std::size_t t = 0; t < matrix_data_ptr_->size(); ++t)
    {
        h_matrix_values_[t] = std::get<2>(matrix_data_ptr_->at(t));
    }
}

void real_dense_solver_cuda::initialize(int matrix_rows, int matrix_columns)
{
    // set the sizes of the system components:
    matrix_cols_ = matrix_columns;
    matrix_rows_ = matrix_rows;

    // clear the containers:
    h_matrix_values_.clear();
    h_rhs_values_.clear();

    // resize the containers to the correct size:
    h_matrix_values_.resize(matrix_cols_ * matrix_rows_);
    h_rhs_values_.resize(matrix_rows_);
}

void real_dense_solver_cuda::solve(std::vector<double> &container, factorization_enum factorization)
{
    populate();

    // get the dimensions:
    std::size_t lda = std::max(matrix_data_ptr_->rows(), matrix_data_ptr_->columns());
    std::size_t m = std::min(matrix_data_ptr_->rows(), matrix_data_ptr_->columns());
    std::size_t ldb = h_rhs_values_.size();

    // step 1: create device containers:
    thrust::device_vector<double> d_matrix_values = h_matrix_values_;
    thrust::device_vector<double> d_rhs_values = h_rhs_values_;
    thrust::device_vector<double> d_solution(m);

    // step 2: cast to raw pointers:
    double *d_mat_vals = thrust::raw_pointer_cast(d_matrix_values.data());
    double *d_rhs_vals = thrust::raw_pointer_cast(d_rhs_values.data());
    double *d_sol = thrust::raw_pointer_cast(d_solution.data());

    // create the helpers:
    real_dense_solver_cuda_helpers helpers;
    helpers.initialize();

    // call the DenseSolverPolicy:
    if (factorization == factorization_enum::QRMethod)
    {
        dense_solver_qr::solve(helpers.get_dense_solver_handle(), helpers.get_cublas_handle(), m, d_mat_vals, lda,
                               d_rhs_vals, d_sol);
    }
    else if (factorization == factorization_enum::LUMethod)
    {
        dense_solver_lu::solve(helpers.get_dense_solver_handle(), helpers.get_cublas_handle(), m, d_mat_vals, lda,
                               d_rhs_vals, d_sol);
    }
    else if (factorization == factorization_enum::CholeskyMethod)
    {
        dense_solver_cholesky::solve(helpers.get_dense_solver_handle(), helpers.get_cublas_handle(), m, d_mat_vals, lda,
                                     d_rhs_vals, d_sol);
    }
    else
    {
        throw std::exception("factorization not known");
    }
    thrust::copy(d_solution.begin(), d_solution.end(), container.begin());
}

container_t const real_dense_solver_cuda::solve(
    factorization_enum factorization) {
    populate();

    // get the dimensions:
    std::size_t lda = std::max(matrix_data_ptr_->rows(), matrix_data_ptr_->columns());
    std::size_t m = std::min(matrix_data_ptr_->rows(), matrix_data_ptr_->columns());
    std::size_t ldb = h_rhs_values_.size();

    // prepare container for solution:
    thrust::host_vector<double> h_solution(m);

    // step 1: create device vectors:
    thrust::device_vector<double> d_matrix_values = h_matrix_values_;
    thrust::device_vector<double> d_rhs_values = h_rhs_values_;
    thrust::device_vector<double> d_solution = h_solution;

    // step 2: cast to raw pointers:
    double *d_mat_vals = thrust::raw_pointer_cast(d_matrix_values.data());
    double *d_rhs_vals = thrust::raw_pointer_cast(d_rhs_values.data());
    double *d_sol = thrust::raw_pointer_cast(d_solution.data());

    // create the helpers:
    real_dense_solver_cuda_helpers helpers;
    helpers.initialize();

    // call the DenseSolverPolicy:
    if (factorization == factorization_enum::QRMethod)
    {
        dense_solver_qr::solve(helpers.get_dense_solver_handle(), helpers.get_cublas_handle(), m, d_mat_vals, lda,
                               d_rhs_vals, d_sol);
    }
    else if (factorization == factorization_enum::LUMethod)
    {
        dense_solver_lu::solve(helpers.get_dense_solver_handle(), helpers.get_cublas_handle(), m, d_mat_vals, lda,
                               d_rhs_vals, d_sol);
    }
    else if (factorization == factorization_enum::CholeskyMethod)
    {
        dense_solver_cholesky::solve(helpers.get_dense_solver_handle(), helpers.get_cublas_handle(), m, d_mat_vals, lda,
                                     d_rhs_vals, d_sol);
    }
    else
    {
        throw std::exception("factorization not known");
    }
    container_t solution(h_solution.size());
    thrust::copy(d_solution.begin(), d_solution.end(), solution.begin());

    return solution;
}

} // namespace lss_dense_solvers
