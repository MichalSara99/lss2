#include "lss_core_cuda_solver.hpp"

#include <cuda_runtime.h>
#include <cusolverSp.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <type_traits>
#include "lss_core_cuda_solver_policy.hpp"

namespace lss_core_cuda_solver
{

using lss_core_cuda_solver_policy::sparse_solver_device;
using lss_core_cuda_solver_policy::sparse_solver_device_cholesky;
using lss_core_cuda_solver_policy::sparse_solver_device_qr;
using lss_core_cuda_solver_policy::sparse_solver_host;
using lss_core_cuda_solver_policy::sparse_solver_host_cholesky;
using lss_core_cuda_solver_policy::sparse_solver_host_lu;
using lss_core_cuda_solver_policy::sparse_solver_host_qr;

real_sparse_solver_cuda<memory_space_enum::Host>::real_sparse_solver_cuda()
{
}


real_sparse_solver_cuda<memory_space_enum::Host>::real_sparse_solver_cuda(
    int system_size)
    : system_size_{system_size}
{
}


real_sparse_solver_cuda<memory_space_enum::Host>::~real_sparse_solver_cuda()
{
}

std::size_t real_sparse_solver_cuda<memory_space_enum::Host>::non_zero_elements() const
{
    return matrix_data_ptr_->size();
}

void real_sparse_solver_cuda<memory_space_enum::Host>::set_rhs(
    container_t const &rhs)
{
    LSS_ASSERT(rhs.size() == system_size_, " rhs has incorrect size");
    copy(h_vector_values_, rhs);
}

void real_sparse_solver_cuda<memory_space_enum::Host>::set_rhs_value(
    std::size_t idx, double value)
{
    LSS_ASSERT((idx < system_size_), "idx is outside range");
    h_vector_values_[idx] = value;
}


void real_sparse_solver_cuda<memory_space_enum::Host>::set_flat_sparse_matrix(flat_matrix matrix)
{
    LSS_ASSERT(matrix.columns() == system_size_, " Incorrect number of columns");
    LSS_ASSERT(matrix.rows() == system_size_, " Incorrect number of rows");
    matrix_data_ptr_ = std::make_unique<flat_matrix>(std::move(matrix));
}


void real_sparse_solver_cuda<lss_enumerations::memory_space_enum::Host>::initialize(std::size_t system_size)
{
    // set the size of the linear system:
    system_size_ = system_size;

    // clear the containers:
    h_matrix_values_.clear();
    h_vector_values_.clear();
    h_column_indices_.clear();
    h_row_counts_.clear();

    // resize the containers to the correct size:
    h_vector_values_.resize(system_size_);
    h_row_counts_.resize((system_size_ + 1));
}


void real_sparse_solver_cuda<lss_enumerations::memory_space_enum::Host>::build_csr()
{
    int const nonZeroSize = non_zero_elements();

    LSS_VERIFY(matrix_data_ptr_, "flat_matrix has not been provided.");
    // CUDA sparse solver is row-major:
    matrix_data_ptr_->sort(flat_matrix_sort_enum::RowMajor);

    h_column_indices_.resize(nonZeroSize);
    h_matrix_values_.resize(nonZeroSize);

    int nElement{0};
    int nRowElement{0};
    int lastRow{0};
    h_row_counts_[nRowElement++] = nElement;

    for (int i = 0; i < nonZeroSize; ++i)
    {
        if (lastRow < std::get<0>(matrix_data_ptr_->at(i)))
        {
            h_row_counts_[nRowElement++] = i;
            lastRow++;
        }

        h_column_indices_[i] = std::get<1>(matrix_data_ptr_->at(i));
        h_matrix_values_[i] = std::get<2>(matrix_data_ptr_->at(i));
    }
    h_row_counts_[nRowElement] = nonZeroSize;
}


void real_sparse_solver_cuda<
    lss_enumerations::memory_space_enum::Host>::solve(container_t &solution, factorization_enum factorization)
{
    build_csr();

    // get the non-zero size:
    int const non_zero_size = non_zero_elements();

    // integer for holding index of row where singularity occurs:
    int singular_idx{0};

    // prepare container for solution:
    thrust::host_vector<double> h_solution(system_size_);

    // get the raw host pointers
    double *h_mat_vals = thrust::raw_pointer_cast(h_matrix_values_.data());
    double *h_rhs_vals = thrust::raw_pointer_cast(h_vector_values_.data());
    double *h_sol = thrust::raw_pointer_cast(h_solution.data());
    int *h_col = thrust::raw_pointer_cast(h_column_indices_.data());
    int *h_row = thrust::raw_pointer_cast(h_row_counts_.data());

    // create the helpers:
    real_sparse_solver_cuda_helpers helpers;
    helpers.initialize();
    // call the sparse_solver_host_policy:
    if (factorization == factorization_enum::QRMethod)
    {
        sparse_solver_host_qr::solve(helpers.get_solver_handle(), helpers.get_matrix_descriptor(), system_size_,
                                            non_zero_size, h_mat_vals, h_row, h_col, h_rhs_vals, 0.0, 0, h_sol,
                                            &singular_idx);
    }
    else if (factorization == factorization_enum::LUMethod)
    {
        sparse_solver_host_lu::solve(helpers.get_solver_handle(), helpers.get_matrix_descriptor(), system_size_,
                                            non_zero_size, h_mat_vals, h_row, h_col, h_rhs_vals, 0.0, 0, h_sol,
                                            &singular_idx);
    }
    else if (factorization == factorization_enum::CholeskyMethod)
    {
        sparse_solver_host_cholesky::solve(helpers.get_solver_handle(), helpers.get_matrix_descriptor(),
                                                  system_size_, non_zero_size, h_mat_vals, h_row, h_col, h_rhs_vals,
                                                  0.0, 0, h_sol, &singular_idx);
    }
    else
    {
        throw std::exception("factorization not known");
    }

    LSS_ASSERT(singular_idx < 0, "Sparse matrix is singular at row: " << singular_idx << "\n");
    copy(solution,h_solution);
}


container_t const real_sparse_solver_cuda<lss_enumerations::memory_space_enum::Host>::solve(factorization_enum
                                                                                                     factorization)
{
    build_csr();

    // get the non-zero size:
    int const non_zero_size = non_zero_elements();

    // integer for holding index of row where singularity occurs:
    int singular_idx{0};

    // prepare container for solution:
    thrust::host_vector<double> h_solution(system_size_);

    // get the raw host pointers
    double *h_mat_vals = thrust::raw_pointer_cast(h_matrix_values_.data());
    double *h_rhs_vals = thrust::raw_pointer_cast(h_vector_values_.data());
    double *h_sol = thrust::raw_pointer_cast(h_solution.data());
    int *h_col = thrust::raw_pointer_cast(h_column_indices_.data());
    int *h_row = thrust::raw_pointer_cast(h_row_counts_.data());

    // create the helpers:
    real_sparse_solver_cuda_helpers helpers;
    helpers.initialize();
    // call the sparse_solver_host_policy:
    if (factorization == factorization_enum::QRMethod)
    {
        sparse_solver_host_qr::solve(helpers.get_solver_handle(), helpers.get_matrix_descriptor(), system_size_,
                                            non_zero_size, h_mat_vals, h_row, h_col, h_rhs_vals, 0.0, 0, h_sol,
                                            &singular_idx);
    }
    else if (factorization == factorization_enum::LUMethod)
    {
        sparse_solver_host_lu::solve(helpers.get_solver_handle(), helpers.get_matrix_descriptor(), system_size_,
                                            non_zero_size, h_mat_vals, h_row, h_col, h_rhs_vals, 0.0, 0, h_sol,
                                            &singular_idx);
    }
    else if (factorization == factorization_enum::CholeskyMethod)
    {
        sparse_solver_host_cholesky::solve(helpers.get_solver_handle(), helpers.get_matrix_descriptor(),
                                                  system_size_, non_zero_size, h_mat_vals, h_row, h_col, h_rhs_vals,
                                                  0.0, 0, h_sol, &singular_idx);
    }
    else
    {
        throw std::exception("factorization not known");
    }

    LSS_ASSERT(singular_idx < 0, "Sparse matrix is singular at row: " << singular_idx << "\n");
    container_t solution(h_solution.size());
    copy(solution,h_solution);
    return solution;
}



// single precision DEVICE specialization

real_sparse_solver_cuda<memory_space_enum::Device>::real_sparse_solver_cuda()
{
}

real_sparse_solver_cuda<memory_space_enum::Device>::real_sparse_solver_cuda(
    int system_size)
    : system_size_{system_size}
{
}

real_sparse_solver_cuda<memory_space_enum::Device>::~real_sparse_solver_cuda()
{
}

std::size_t real_sparse_solver_cuda<memory_space_enum::Device>::non_zero_elements() const
{
    LSS_VERIFY(matrix_data_ptr_, "flat_matrix has not been provided.");
    return matrix_data_ptr_->size();
}

void real_sparse_solver_cuda<memory_space_enum::Device>::set_rhs(
    container_t const &rhs)
{
    LSS_ASSERT(rhs.size() == system_size_, " rhs has incorrect size");
    copy(h_vector_values_,rhs);
}


void real_sparse_solver_cuda<memory_space_enum::Device>::set_rhs_value(
    std::size_t idx, double value)
{
    LSS_ASSERT((idx < system_size_), "idx is outside range");
    h_vector_values_[idx] = value;
}


void real_sparse_solver_cuda<memory_space_enum::Device>::set_flat_sparse_matrix(flat_matrix matrix)
{
    LSS_ASSERT(matrix.columns() == system_size_, " Incorrect number of columns");
    LSS_ASSERT(matrix.rows() == system_size_, " Incorrect number of rows");
    matrix_data_ptr_ = std::make_unique<flat_matrix>(std::move(matrix));
}


void real_sparse_solver_cuda<lss_enumerations::memory_space_enum::Device>::initialize(std::size_t system_size)
{
    // set the size of the linear system:
    system_size_ = system_size;

    // clear the containers:
    h_matrix_values_.clear();
    h_vector_values_.clear();
    h_column_indices_.clear();
    h_row_counts_.clear();

    // resize the containers to the correct size:
    h_vector_values_.resize(system_size_);
    h_row_counts_.resize((system_size_ + 1));
}


void real_sparse_solver_cuda<lss_enumerations::memory_space_enum::Device>::build_csr()
{
    int const non_zero_size = non_zero_elements();

    LSS_VERIFY(matrix_data_ptr_, "flat_matrix has not been provided.");
    // CUDA sparse solver is row-major:
    matrix_data_ptr_->sort(flat_matrix_sort_enum::RowMajor);

    h_column_indices_.resize(non_zero_size);
    h_matrix_values_.resize(non_zero_size);

    int nElement{0};
    int nRowElement{0};
    int lastRow{0};
    h_row_counts_[nRowElement++] = nElement;

    for (std::size_t i = 0; i < non_zero_size; ++i)
    {
        if (lastRow < std::get<0>(matrix_data_ptr_->at(i)))
        {
            h_row_counts_[nRowElement++] = i;
            lastRow++;
        }

        h_column_indices_[i] = std::get<1>(matrix_data_ptr_->at(i));
        h_matrix_values_[i] = std::get<2>(matrix_data_ptr_->at(i));
    }
    h_row_counts_[nRowElement] = non_zero_size;
}


void real_sparse_solver_cuda<
    lss_enumerations::memory_space_enum::Device>::solve(container_t &solution, factorization_enum factorization)
{
    build_csr();

    // get the non-zero size:
    int const non_zero_size = non_zero_elements();

    // integer for holding index of row where singularity occurs:
    int singular_idx{0};

    // prepare container for solution:
    thrust::device_vector<double> d_solution(system_size_);

    // copy to the device constainers:
    thrust::device_vector<double> d_matrix_values = h_matrix_values_;
    thrust::device_vector<double> d_vector_values = h_vector_values_;
    thrust::device_vector<int> d_column_indices = h_column_indices_;
    thrust::device_vector<int> d_row_counts = h_row_counts_;

    // get the raw host pointers
    double *d_mat_vals = thrust::raw_pointer_cast(d_matrix_values.data());
    double *d_rhs_vals = thrust::raw_pointer_cast(d_vector_values.data());
    double *d_sol = thrust::raw_pointer_cast(d_solution.data());
    int *d_col = thrust::raw_pointer_cast(d_column_indices.data());
    int *d_row = thrust::raw_pointer_cast(d_row_counts.data());

    // create the helpers:
    real_sparse_solver_cuda_helpers helpers;
    helpers.initialize();
    // call the sparse_solver_device_policy:
    if (factorization == factorization_enum::QRMethod)
    {
        sparse_solver_device_qr::solve(helpers.get_solver_handle(), helpers.get_matrix_descriptor(),
                                              system_size_, non_zero_size, d_mat_vals, d_row, d_col, d_rhs_vals, 0.0, 0,
                                              d_sol, &singular_idx);
    }
    else if (factorization == factorization_enum::CholeskyMethod)
    {
        sparse_solver_device_cholesky::solve(helpers.get_solver_handle(), helpers.get_matrix_descriptor(),
                                                    system_size_, non_zero_size, d_mat_vals, d_row, d_col, d_rhs_vals,
                                                    0.0, 0, d_sol, &singular_idx);
    }
    else if (factorization == factorization_enum::LUMethod)
    {
        throw std::exception("factorization not supported on device.");
    }
    else
    {
        throw std::exception("factorization not known");
    }

    LSS_ASSERT(singular_idx < 0, "Sparse matrix is singular at row: " << singular_idx << "\n");
    copy(solution,d_solution);
}


container_t const real_sparse_solver_cuda<
    lss_enumerations::memory_space_enum::Device>::solve(factorization_enum  factorization)
{
    build_csr();

    // get the non-zero size:
    int const non_zero_size = non_zero_elements();

    // integer for holding index of row where singularity occurs:
    int singular_idx{0};

    // prepare container for solution:
    thrust::device_vector<double> d_solution(system_size_);

    // copy to the device constainers:
    thrust::device_vector<double> d_matrix_values = h_matrix_values_;
    thrust::device_vector<double> d_vector_values = h_vector_values_;
    thrust::device_vector<int> d_column_indices = h_column_indices_;
    thrust::device_vector<int> d_row_counts = h_row_counts_;

    // get the raw host pointers
    double *d_mat_vals = thrust::raw_pointer_cast(d_matrix_values.data());
    double *d_rhs_vals = thrust::raw_pointer_cast(d_vector_values.data());
    double *d_sol = thrust::raw_pointer_cast(d_solution.data());
    int *d_col = thrust::raw_pointer_cast(d_column_indices.data());
    int *d_row = thrust::raw_pointer_cast(d_row_counts.data());

    // create the helpers:
    real_sparse_solver_cuda_helpers helpers;
    helpers.initialize();
    // call the sparse_solver_device_policy:
    if (factorization == factorization_enum::QRMethod)
    {
        sparse_solver_device_qr::solve(helpers.get_solver_handle(), helpers.get_matrix_descriptor(),
                                              system_size_, non_zero_size, d_mat_vals, d_row, d_col, d_rhs_vals, 0.0, 0,
                                              d_sol, &singular_idx);
    }
    else if (factorization == factorization_enum::CholeskyMethod)
    {
        sparse_solver_device_cholesky::solve(helpers.get_solver_handle(), helpers.get_matrix_descriptor(),
                                                    system_size_, non_zero_size, d_mat_vals, d_row, d_col, d_rhs_vals,
                                                    0.0, 0, d_sol, &singular_idx);
    }
    else if (factorization == factorization_enum::LUMethod)
    {
        throw std::exception("factorization not supported on device.");
    }
    else
    {
        throw std::exception("factorization not known");
    }

    LSS_ASSERT(singular_idx < 0, "Sparse matrix is singular at row: " << singular_idx << "\n");
    container_t solution(system_size_);
    copy(solution,d_solution);
    return solution;
}

} // namespace lss_core_cuda_solver
