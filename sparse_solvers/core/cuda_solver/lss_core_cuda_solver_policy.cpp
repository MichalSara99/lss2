#include "lss_core_cuda_solver_policy.hpp"

namespace lss_core_cuda_solver_policy
{

void sparse_solver_host_qr::solve(cusolverSpHandle_t handle, cusparseMatDescr_t mat_desc, int const system_size,
                                  int const non_zero_size, double const *h_mat_vals, int const *h_row_counts,
                                  int const *h_col_indices, double const *h_rhs, double tol, int reorder,
                                  double *h_solution, int *singularity)
{
    CUSOLVER_STATUS(cusolverSpDcsrlsvqrHost(handle, system_size, non_zero_size, mat_desc, h_mat_vals, h_row_counts,
                                            h_col_indices, h_rhs, tol, reorder, h_solution, singularity));
}

void sparse_solver_host_lu::solve(cusolverSpHandle_t handle, cusparseMatDescr_t mat_desc, int const system_size,
                                  int const non_zero_size, double const *h_mat_vals, int const *h_row_counts,
                                  int const *h_col_indices, double const *h_rhs, double tol, int reorder,
                                  double *h_solution, int *singularity)
{
    CUSOLVER_STATUS(cusolverSpDcsrlsvluHost(handle, system_size, non_zero_size, mat_desc, h_mat_vals, h_row_counts,
                                            h_col_indices, h_rhs, tol, reorder, h_solution, singularity));
}

void sparse_solver_host_cholesky::solve(cusolverSpHandle_t handle, cusparseMatDescr_t mat_desc, int const system_size,
                                        int const non_zero_size, double const *h_mat_vals, int const *h_row_counts,
                                        int const *h_col_indices, double const *h_rhs, double tol, int reorder,
                                        double *h_solution, int *singularity)
{
    CUSOLVER_STATUS(cusolverSpDcsrlsvcholHost(handle, system_size, non_zero_size, mat_desc, h_mat_vals, h_row_counts,
                                              h_col_indices, h_rhs, tol, reorder, h_solution, singularity));
}

void sparse_solver_device_qr::solve(cusolverSpHandle_t handle, cusparseMatDescr_t mat_desc, int const system_size,
                                    int const non_zero_size, double const *d_mat_vals, int const *d_row_counts,
                                    int const *d_col_indices, double const *d_rhs, double tol, int reorder,
                                    double *d_solution, int *singularity)
{
    CUSOLVER_STATUS(cusolverSpDcsrlsvqr(handle, system_size, non_zero_size, mat_desc, d_mat_vals, d_row_counts,
                                        d_col_indices, d_rhs, tol, reorder, d_solution, singularity));
}

void sparse_solver_device_cholesky::solve(cusolverSpHandle_t handle, cusparseMatDescr_t mat_desc, int const system_size,
                                          int const non_zero_size, double const *d_mat_vals, int const *d_row_counts,
                                          int const *d_col_indices, double const *d_rhs, double tol, int reorder,
                                          double *d_solution, int *singularity)
{
    CUSOLVER_STATUS(cusolverSpDcsrlsvchol(handle, system_size, non_zero_size, mat_desc, d_mat_vals, d_row_counts,
                                          d_col_indices, d_rhs, tol, reorder, d_solution, singularity));
}
} // namespace lss_core_cuda_solver_policy
