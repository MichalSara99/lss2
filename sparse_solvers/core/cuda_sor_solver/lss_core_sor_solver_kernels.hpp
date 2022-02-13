#pragma once
#if !defined(_LSS_CORE_SOR_SOLVER_KERNELS_HPP_)
#define _LSS_CORE_SOR_SOLVER_KERNELS_HPP_

#include "../../../common/lss_utility.hpp"
#include <device_launch_parameters.h>

#define THREADS_PER_BLOCK 256

namespace lss_core_sor_solver
{

__global__ void sor_kernel(double const *matrix_vals, std::size_t row_size, std::size_t const *column_idx,
                           double const *rhv_vals, double const *diagonal_vals, std::size_t const *row_start_idx,
                           std::size_t const *row_end_idx, double const omega, double *sol, double *new_sol,
                           double *errors)
{
    std::size_t const tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid >= row_size)
        return;

    double sigma_1{0.0};
    double sigma_2{0.0};
    double mat_val{0.0};
    const double one = static_cast<double>(1.0);
    double const diag = diagonal_vals[tid];
    std::size_t col_idx{};
    std::size_t const start_idx = row_start_idx[tid];
    std::size_t const end_idx = row_end_idx[tid];
    for (std::size_t c = start_idx; c <= end_idx; ++c)
    {
        mat_val = matrix_vals[c];
        col_idx = column_idx[c];
        if (col_idx < tid)
            sigma_1 += mat_val * new_sol[col_idx];
        if (col_idx > tid)
            sigma_2 += mat_val * sol[col_idx];
    }
    new_sol[tid] = (one - omega) * sol[tid] + ((omega / diag) * (rhv_vals[tid] - sigma_1 - sigma_2));
    errors[tid] = (new_sol[tid] - sol[tid]) * (new_sol[tid] - sol[tid]);
}

} // namespace lss_core_sor_solver

#endif ///_LSS_CORE_SOR_SOLVER_KERNELS_HPP_
