#if !defined(_LSS_SOR_SOLVER_CORE_HPP_)
#define _LSS_SOR_SOLVER_CORE_HPP_

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include <numeric>

#include "../../../common/lss_utility.hpp"

namespace lss_core_sor_solver
{

using lss_utility::NaN;
using lss_utility::sptr_t;
using lss_utility::uptr_t;

class sor_solver_core
{
  private:
    thrust::device_vector<double> matrix_vals_;
    thrust::device_vector<std::size_t> non_zero_column_idx_;
    thrust::device_vector<double> rhs_vals_;
    thrust::device_vector<double> diagonal_vals_;
    thrust::device_vector<std::size_t> row_start_idx_;
    thrust::device_vector<std::size_t> row_end_idx_;
    std::size_t nrows_;

    explicit sor_solver_core();

  public:
    explicit sor_solver_core(thrust::device_vector<double> const &matrix_values,
                             thrust::device_vector<std::size_t> const &non_zero_column_idx, std::size_t const &nrows,
                             thrust::device_vector<double> const &rhs_values,
                             thrust::device_vector<double> const &diagonal_values,
                             thrust::device_vector<std::size_t> const &row_start_idx,
                             thrust::device_vector<std::size_t> const &row_end_idx);

    ~sor_solver_core();

    sor_solver_core(sor_solver_core const &) = delete;
    sor_solver_core(sor_solver_core &&) = delete;
    sor_solver_core &operator=(sor_solver_core const &) = delete;
    sor_solver_core &operator=(sor_solver_core &&) = delete;

    void operator()(thrust::device_vector<double> &solution, thrust::device_vector<double> &next_solution,
                    thrust::device_vector<double> &errors, double omega);
};
} // namespace lss_core_sor_solver

#endif ///_LSS_SOR_SOLVER_CORE_HPP_
