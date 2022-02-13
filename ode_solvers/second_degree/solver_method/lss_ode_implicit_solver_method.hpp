/**

    @file      lss_ode_implicit_solver_method.hpp
    @brief     ODE implicit solver method
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_GENERAL_ODE_IMPLICIT_SOLVER_METHOD_HPP_)
#define _LSS_GENERAL_ODE_IMPLICIT_SOLVER_METHOD_HPP_

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../discretization/lss_discretization.hpp"
#include "../../../discretization/lss_grid.hpp"
#include "../../../discretization/lss_grid_config.hpp"
#include "../../../ode_solvers/second_degree/implicit_coefficients/lss_ode_implicit_coefficients.hpp"
#include "../../../sparse_solvers/tridiagonal/lss_tridiagonal_solver.hpp"

namespace lss_ode_solvers
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_1d_ptr;
using lss_tridiagonal_solver::tridiagonal_solver_ptr;
using lss_utility::container_t;

/**
wave_implicit_solver_method object
*/
class ode_implicit_solver_method
{

  private:
    // solvers:
    tridiagonal_solver_ptr solveru_ptr_;
    // scheme coefficients:
    ode_implicit_coefficients_ptr coefficients_;
    grid_config_1d_ptr grid_cfg_;
    // container:
    container_t low_, diag_, high_, rhs_;

    explicit ode_implicit_solver_method() = delete;

    void initialize();

    void split(container_t &low, container_t &diag, container_t &high);

  public:
    explicit ode_implicit_solver_method(tridiagonal_solver_ptr const &solver_ptr,
                                        ode_implicit_coefficients_ptr const &coefficients,
                                        grid_config_1d_ptr const &grid_config);

    ~ode_implicit_solver_method();

    ode_implicit_solver_method(ode_implicit_solver_method const &) = delete;
    ode_implicit_solver_method(ode_implicit_solver_method &&) = delete;
    ode_implicit_solver_method &operator=(ode_implicit_solver_method const &) = delete;
    ode_implicit_solver_method &operator=(ode_implicit_solver_method &&) = delete;

    void solve(boundary_1d_pair const &boundary_pair, container_t &solution);

    void solve(boundary_1d_pair const &boundary_pair, std::function<double(double)> const &source,
               container_t &solution);
};

} // namespace lss_ode_solvers

#endif ///_LSS_GENERAL_ODE_IMPLICIT_SOLVER_METHOD_HPP_
