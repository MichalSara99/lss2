/**

    @file      lss_heat_implicit_solver_method.hpp
    @brief     Implicit Solver method
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_IMPLICIT_SOLVER_METHOD_HPP_)
#define _LSS_HEAT_IMPLICIT_SOLVER_METHOD_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../../../../sparse_solvers/tridiagonal/lss_tridiagonal_solver.hpp"
#include "../implicit_coefficients/lss_heat_coefficients.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_grids::grid_config_1d_ptr;
using lss_utility::container_t;
using lss_utility::sptr_t;

class implicit_heat_scheme
{

  public:
    static void rhs(heat_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg, container_t const &input,
                    boundary_1d_pair const &boundary_pair, double const &time, container_t &solution);

    static void rhs_source(heat_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                           container_t const &input, container_t const &inhom_input,
                           container_t const &inhom_input_next, boundary_1d_pair const &boundary_pair,
                           double const &time, container_t &solution);
};

/**
heat_implicit_solver_method object
*/
class heat_implicit_solver_method
{
  private:
    // private constats:
    const double cone_ = 1.0;
    const double ctwo_ = 2.0;
    // solvers:
    lss_tridiagonal_solver::tridiagonal_solver_ptr solveru_ptr_;
    // scheme coefficients:
    heat_coefficients_ptr coefficients_;
    grid_config_1d_ptr grid_cfg_;
    // prepare containers:
    container_t low_, diag_, high_;
    container_t source_, source_next_;
    container_t rhs_;

    explicit heat_implicit_solver_method() = delete;

    void initialize(bool is_heat_sourse_set);

    void split(double const &time, container_t &low, container_t &diag, container_t &high);

  public:
    explicit heat_implicit_solver_method(lss_tridiagonal_solver::tridiagonal_solver_ptr const &solver_ptr,
                                         heat_coefficients_ptr const &coefficients,
                                         grid_config_1d_ptr const &grid_config, bool is_heat_sourse_set);

    ~heat_implicit_solver_method();

    heat_implicit_solver_method(heat_implicit_solver_method const &) = delete;
    heat_implicit_solver_method(heat_implicit_solver_method &&) = delete;
    heat_implicit_solver_method &operator=(heat_implicit_solver_method const &) = delete;
    heat_implicit_solver_method &operator=(heat_implicit_solver_method &&) = delete;

    void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, double const &time,
               container_t &solution);

    void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, double const &time,
               double const &next_time, std::function<double(double, double)> const &heat_source,
               container_t &solution);
};

using heat_implicit_solver_method_ptr = sptr_t<heat_implicit_solver_method>;

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_IMPLICIT_SOLVER_METHOD_HPP_
