/**

    @file      lss_double_sweep_solver.hpp
    @brief     Tridiagonal Double Sweep Solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_DOUBLE_SWEEP_SOLVER_HPP_)
#define _LSS_DOUBLE_SWEEP_SOLVER_HPP_

#pragma warning(disable : 4244)

#include <tuple>
#include <type_traits>
#include <vector>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "../lss_tridiagonal_solver.hpp"
#include "lss_double_sweep_boundary.hpp"

namespace lss_double_sweep_solver
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_3d_pair;
using lss_enumerations::factorization_enum;
using lss_utility::container_t;
using lss_utility::sptr_t;

/**
 * fdm_double_sweep_solver
 */

class double_sweep_solver : public lss_tridiagonal_solver::tridiagonal_solver
{
  private:
    container_t L_, K_;
    double_sweep_boundary_ptr dss_boundary_;

    void kernel(boundary_1d_pair const &boundary, container_t &solution, factorization_enum factorization,
                double time) override;

    void kernel(boundary_2d_pair const &boundary, container_t &solution, factorization_enum factorization, double time,
                double space_arg) override;

    void kernel(boundary_3d_pair const &boundary, container_t &solution, factorization_enum factorization, double time,
                double space_1_arg, double space_2_arg) override;

    void initialize();

    explicit double_sweep_solver() = delete;

  public:
    explicit double_sweep_solver(std::size_t discretization_size);

    ~double_sweep_solver();
};

using double_sweep_solver_ptr = sptr_t<double_sweep_solver>;

} // namespace lss_double_sweep_solver

#endif ///_LSS_DOUBLE_SWEEP_SOLVER_HPP_
