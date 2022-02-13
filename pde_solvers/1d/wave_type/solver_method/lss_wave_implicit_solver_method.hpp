/**

    @file      lss_wave_implicit_solver_method.hpp
    @brief     Implicit solver method for wave problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_WAVE_IMPLICIT_SOLVER_METHOD_HPP_)
#define _LSS_WAVE_IMPLICIT_SOLVER_METHOD_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../../../../sparse_solvers/tridiagonal/lss_tridiagonal_solver.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../implicit_coefficients/lss_wave_implicit_coefficients.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_grids::grid_config_1d;
using lss_utility::container_t;

class implicit_wave_scheme
{

  public:
    static void rhs(wave_implicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                    container_t const &input_0, container_t const &input_1, boundary_1d_pair const &boundary_pair,
                    double const &time, container_t &solution);

    static void rhs_source(wave_implicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                           container_t const &input_0, container_t const &input_1, container_t const &inhom_input,
                           boundary_1d_pair const &boundary_pair, double const &time, container_t &solution);

    static void rhs_initial(wave_implicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                            container_t const &input_0, container_t const &input_1,
                            boundary_1d_pair const &boundary_pair, double const &time, container_t &solution);

    static void rhs_initial_source(wave_implicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                                   container_t const &input_0, container_t const &input_1,
                                   container_t const &inhom_input, boundary_1d_pair const &boundary_pair,
                                   double const &time, container_t &solution);

    static void rhs_terminal(wave_implicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                             container_t const &input_0, container_t const &input_1,
                             boundary_1d_pair const &boundary_pair, double const &time, container_t &solution);

    static void rhs_terminal_source(wave_implicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                                    container_t const &input_0, container_t const &input_1,
                                    container_t const &inhom_input, boundary_1d_pair const &boundary_pair,
                                    double const &time, container_t &solution);
};

/**
wave_implicit_solver_method object
*/
class wave_implicit_solver_method
{

  private:
    // constant coeffs:
    const double ctwo_ = 2.0;
    const double cone_ = 1.0;
    // solvers:
    lss_tridiagonal_solver::tridiagonal_solver_ptr solveru_ptr_;
    // scheme coefficients:
    wave_implicit_coefficients_ptr coefficients_;
    grid_config_1d_ptr grid_cfg_;
    // container:
    container_t low_, diag_, high_;
    container_t rhs_, source_;

    explicit wave_implicit_solver_method() = delete;

    void initialize(bool is_wave_source_set);

    void split_0(double time, container_t &low_0, container_t &diag_0, container_t &high_0);

    void split_1(double time, container_t &low_1, container_t &diag_1, container_t &high_1);

  public:
    explicit wave_implicit_solver_method(lss_tridiagonal_solver::tridiagonal_solver_ptr const &solver_ptr,
                                         wave_implicit_coefficients_ptr const &coefficients,
                                         grid_config_1d_ptr const &grid_config, bool is_wave_source_set);

    ~wave_implicit_solver_method();

    wave_implicit_solver_method(wave_implicit_solver_method const &) = delete;
    wave_implicit_solver_method(wave_implicit_solver_method &&) = delete;
    wave_implicit_solver_method &operator=(wave_implicit_solver_method const &) = delete;
    wave_implicit_solver_method &operator=(wave_implicit_solver_method &&) = delete;

    void solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                       boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                       container_t &solution);

    void solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                       boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                       std::function<double(double, double)> const &wave_source, container_t &solution);

    void solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                        boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                        container_t &solution);

    void solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                        boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                        std::function<double(double, double)> const &wave_source, container_t &solution);

    void solve(container_t &prev_solution_0, container_t &prev_solution_1, boundary_1d_pair const &boundary_pair,
               double const &time, double const &next_time, container_t &solution);

    void solve(container_t &prev_solution_0, container_t &prev_solution_1, boundary_1d_pair const &boundary_pair,
               double const &time, double const &next_time, std::function<double(double, double)> const &wave_source,
               container_t &solution);
};

using wave_implicit_solver_method_ptr = sptr_t<wave_implicit_solver_method>;

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_WAVE_IMPLICIT_SOLVER_METHOD_HPP_
