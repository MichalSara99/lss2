#include "lss_core_sor_solver_kernels.hpp"
#include "lss_sor_solver_core.hpp"

namespace lss_core_sor_solver {


    sor_solver_core::sor_solver_core() {}

    sor_solver_core::sor_solver_core(
        thrust::device_vector<double> const &matrix_values,
        thrust::device_vector<std::size_t> const &non_zero_column_idx,
        std::size_t const &nrows,
        thrust::device_vector<double> const &rhs_values,
        thrust::device_vector<double> const &diagonal_values,
        thrust::device_vector<std::size_t> const &row_start_idx,
        thrust::device_vector<std::size_t> const &row_end_idx)
        : matrix_vals_{matrix_values},
          non_zero_column_idx_{non_zero_column_idx},
          nrows_{nrows},
          rhs_vals_{rhs_values},
          diagonal_vals_{diagonal_values},
          row_start_idx_{row_start_idx},
          row_end_idx_{row_end_idx} {}


    sor_solver_core::~sor_solver_core() {}

    void sor_solver_core::operator()(
        thrust::device_vector<double> &solution,
        thrust::device_vector<double> &next_solution,
        thrust::device_vector<double> &errors, double omega) {
      unsigned long long const threads_per_block = THREADS_PER_BLOCK;
      unsigned long long const blocks_per_grid =
          (nrows_ + threads_per_block - 1) / threads_per_block;

      double const *mat_data_ptr = thrust::raw_pointer_cast(&matrix_vals_[0]);
      std::size_t const *non_zero_col_idx_data_ptr =
          thrust::raw_pointer_cast(&non_zero_column_idx_[0]);
      double const *rhs_data_ptr = thrust::raw_pointer_cast(&rhs_vals_[0]);
      double const *diag_data_ptr =
          thrust::raw_pointer_cast(&diagonal_vals_[0]);
      std::size_t const *row_start_idx_data_ptr =
          thrust::raw_pointer_cast(&row_start_idx_[0]);
      std::size_t const *row_end_idx_data_ptr =
          thrust::raw_pointer_cast(&row_end_idx_[0]);
      double *sol_data_ptr = thrust::raw_pointer_cast(&solution[0]);
      double *next_sol_data_ptr = thrust::raw_pointer_cast(&next_solution[0]);
      double *errors_data_ptr = thrust::raw_pointer_cast(&errors[0]);

      sor_kernel<<<threads_per_block, blocks_per_grid>>>(
              mat_data_ptr, nrows_, non_zero_col_idx_data_ptr, rhs_data_ptr,
              diag_data_ptr, row_start_idx_data_ptr, row_end_idx_data_ptr,
              omega, sol_data_ptr, next_sol_data_ptr, errors_data_ptr);
    }


}




