/**

    @file      lss_wave_euler_solver_method.hpp
    @brief     Euler solver method for wave problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_WAVE_EULER_SOLVER_METHOD_HPP_)
#define _LSS_WAVE_EULER_SOLVER_METHOD_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../explicit_coefficients/lss_wave_explicit_coefficients.hpp"
#include "lss_wave_explicit_solver_method.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_grids::grid_config_1d_ptr;
using lss_utility::container_t;

/**
    explicit_wave_scheme object
 */
class explicit_wave_scheme
{

  public:
    static void rhs(wave_explicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                    container_t const &input_0, container_t const &input_1, boundary_1d_pair const &boundary_pair,
                    double const &time, container_t &solution);

    static void rhs_source(wave_explicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                           container_t const &input_0, container_t const &input_1, container_t const &inhom_input,
                           boundary_1d_pair const &boundary_pair, double const &time, container_t &solution);

    static void rhs_initial(wave_explicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                            container_t const &input_0, container_t const &input_1,
                            boundary_1d_pair const &boundary_pair, double const &time, container_t &solution);

    static void rhs_initial_source(wave_explicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                                   container_t const &input_0, container_t const &input_1,
                                   container_t const &inhom_input, boundary_1d_pair const &boundary_pair,
                                   double const &time, container_t &solution);

    static void rhs_terminal(wave_explicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                             container_t const &input_0, container_t const &input_1,
                             boundary_1d_pair const &boundary_pair, double const &time, container_t &solution);

    static void rhs_terminal_source(wave_explicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                                    container_t const &input_0, container_t const &input_1,
                                    container_t const &inhom_input, boundary_1d_pair const &boundary_pair,
                                    double const &time, container_t &solution);
};

/**
wave_euler_solver_method object
*/
class wave_euler_solver_method : public wave_explicit_solver_method
{

  private:
    // scheme coefficients:
    wave_explicit_coefficients_ptr coefficients_;
    // containers:
    container_t source_;

    explicit wave_euler_solver_method() = delete;

    void initialize(bool is_wave_source_set);

  public:
    explicit wave_euler_solver_method(wave_explicit_coefficients_ptr const &coefficients,
                                      grid_config_1d_ptr const &grid_config, bool is_wave_source_set);

    ~wave_euler_solver_method();

    wave_euler_solver_method(wave_euler_solver_method const &) = delete;
    wave_euler_solver_method(wave_euler_solver_method &&) = delete;
    wave_euler_solver_method &operator=(wave_euler_solver_method const &) = delete;
    wave_euler_solver_method &operator=(wave_euler_solver_method &&) = delete;

    void solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                       boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                       container_t &solution) override;

    void solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                       boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                       std::function<double(double, double)> const &wave_source, container_t &solution) override;

    void solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                        boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                        container_t &solution) override;

    void solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                        boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                        std::function<double(double, double)> const &wave_source, container_t &solution) override;

    void solve(container_t &prev_solution_0, container_t &prev_solution_1, boundary_1d_pair const &boundary_pair,
               double const &time, double const &next_time, container_t &solution) override;

    void solve(container_t &prev_solution_0, container_t &prev_solution_1, boundary_1d_pair const &boundary_pair,
               double const &time, double const &next_time, std::function<double(double, double)> const &wave_source,
               container_t &solution) override;
};

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_WAVE_EULER_SOLVER_METHOD_HPP_
